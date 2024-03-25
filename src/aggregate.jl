# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Aggregate(domain, var₁ => agg₁, var₂ => agg₂, ..., varₙ => aggₙ)
    Aggregate([g₁, g₂, ..., gₙ], var₁ => agg₁, var₂ => agg₂, ..., varₙ => aggₙ)

Aggregate variables `var₁`, `var₂`, ..., `varₙ` over geospatial `domain` using
aggregation functions `agg₁`, `agg₂`, ..., `aggₙ`. Alternatively, aggregate
variables over geometries `g₁`, `g₂`, ..., `gₙ`. Default aggregation function
is `mean` for continuous variables and `first` otherwise.

# Examples

```julia
Aggregate(domain, 1 => last, 2 => maximum)
Aggregate(domain, :a => first, :b => minimum)
Aggregate(domain, "a" => last, "b" => maximum)
Aggregate(geoms, "a" => last, "b" => maximum)
```
"""
struct Aggregate{D<:Domain,S<:ColumnSelector} <: TableTransform
  domain::D
  selector::S
  aggfuns::Vector{Function}
end

Aggregate(domain::Domain) = Aggregate(domain, NoneSelector(), Function[])
Aggregate(domain::Domain, pairs::Pair{C,<:Function}...) where {C<:Column} =
  Aggregate(domain, selector(first.(pairs)), collect(Function, last.(pairs)))
Aggregate(geoms::AbstractVector{<:Geometry}, args...) = Aggregate(GeometrySet(geoms), args...)

isrevertible(::Type{<:Aggregate}) = false

function apply(transform::Aggregate, geotable::AbstractGeoTable)
  gtb = _adjustunits(geotable)
  dom = domain(gtb)
  tab = values(gtb)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  # aggregation functions
  svars = transform.selector(vars)
  agg = Dict(zip(svars, transform.aggfuns))
  for var in vars
    if !haskey(agg, var)
      v = Tables.getcolumn(cols, var)
      agg[var] = _defaultagg(v)
    end
  end

  # find target indices of each row
  tdom = transform.domain
  knn = KNearestSearch(tdom, 1)
  inds = _tmap(1:nelements(dom)) do i
    first(search(centroid(dom, i), knn))
  end

  # group rows with the same target indices
  tinds = 1:nelements(tdom)
  group = Dict(tind => Int[] for tind in tinds)
  for (i, tind) in enumerate(inds)
    push!(group[tind], i)
  end

  # perform aggregation with repeated indices
  function aggvar(var)
    v = Tables.getcolumn(cols, var)
    map(tinds) do tind
      sinds = group[tind]
      agg[var](v[sinds])
    end
  end

  # construct new table
  𝒯 = (; (var => aggvar(var) for var in vars)...)
  newtab = 𝒯 |> Tables.materializer(tab)

  # new spatial data
  newgeotable = georef(newtab, tdom)

  newgeotable, nothing
end
