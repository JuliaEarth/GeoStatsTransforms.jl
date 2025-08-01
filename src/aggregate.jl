# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Aggregate(domain, var₁ => agg₁, var₂ => agg₂, ..., varₙ => aggₙ)

Aggregate variables `var₁`, `var₂`, ..., `varₙ` over geospatial `domain` using
aggregation functions `agg₁`, `agg₂`, ..., `aggₙ`.

    Aggregate([g₁, g₂, ..., gₙ], var₁ => agg₁, var₂ => agg₂, ..., varₙ => aggₙ)

Alternatively, aggregate variables over geometries `g₁`, `g₂`, ..., `gₙ`.

Default aggregation function is `mean` for continuous variables and `first` otherwise.

## Examples

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
  gtb = geotable |> AbsoluteUnits()
  table = values(gtb)
  cols = Tables.columns(table)
  vars = Tables.columnnames(cols)

  # aggregation functions
  svars = transform.selector(vars)
  aggfun = Dict(zip(svars, transform.aggfuns))
  for var in vars
    if !haskey(aggfun, var)
      vals = Tables.getcolumn(cols, var)
      aggfun[var] = _defaultagg(vals)
    end
  end

  # source and target domains
  sdom = domain(gtb)
  tdom = transform.domain

  # perform aggregation
  newcols = _aggregate(sdom, tdom, cols, vars, aggfun)
  newtable = (; newcols...) |> Tables.materializer(table)

  georef(newtable, tdom), nothing
end

function _aggregate(sdom, tdom, cols, vars, aggfun)
  if sdom isa Grid && tdom isa Grid && extrema(sdom) == extrema(tdom) && all(iszero, size(sdom) .% size(tdom))
    # we have two grids overlaid, and can rely on
    # tiled iteration for efficient aggregation
    _gridagg(sdom, tdom, cols, vars, aggfun)
  else
    # general case with knn search
    _knnagg(sdom, tdom, cols, vars, aggfun)
  end
end

function _gridagg(sdom, tdom, cols, vars, aggfun)
  # determine tile size for tiled iteration
  tilesize = size(sdom) .÷ size(tdom)
  if any(<(1), tilesize)
    throw(ArgumentError("cannot aggregate a coarse grid over a fine grid"))
  end

  # perform aggregation
  map(vars) do var
    svals = Tables.getcolumn(cols, var)
    array = reshape(svals, size(sdom))
    titer = TileIterator(axes(array), tilesize)
    tvals = tmap(titer) do sinds
      aggfun[var](array[sinds...])
    end |> vec
    var => tvals
  end
end

function _knnagg(sdom, tdom, cols, vars, aggfun)
  # find nearest elements in target domain
  knn = KNearestSearch(tdom, 1)
  near = tmap(1:nelements(sdom)) do i
    first(search(centroid(sdom, i), knn))
  end

  # map target element to source elements
  group = Dict(tind => Int[] for tind in 1:nelements(tdom))
  for (sind, tind) in enumerate(near)
    push!(group[tind], sind)
  end

  # perform aggregation
  map(vars) do var
    svals = Tables.getcolumn(cols, var)
    tvals = tmap(1:nelements(tdom)) do tind
      aggfun[var](svals[group[tind]])
    end
    var => tvals
  end
end
