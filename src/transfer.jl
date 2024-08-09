# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Transfer(domain)
    Transfer([g₁, g₂, ..., gₙ])

Transfer variables `var₁`, `var₂`, ..., `varₙ` from source domain to target `domain`.
Alternatively, transfer variables to geometries `g₁`, `g₂`, ..., `gₙ`.

# Examples

```julia
Transfer(CartesianGrid(10, 10))
Transfer(rand(Point, 100))
```
"""
struct Transfer{D<:Domain} <: TableTransform
  domain::D
end

Transfer(geoms::AbstractVector{<:Geometry}) = Transfer(GeometrySet(geoms))

isrevertible(::Type{<:Transfer}) = false

function apply(transform::Transfer, geotable::AbstractGeoTable)
  gtb = _adjustunits(geotable)
  dom = domain(gtb)
  tab = values(gtb)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  # find target indices of each row
  tdom = transform.domain
  knn = KNearestSearch(dom, 1)
  inds = tmap(1:nelements(tdom)) do i
    first(search(centroid(tdom, i), knn))
  end

  # generate variables for target domain
  function genvar(var)
    v = Tables.getcolumn(cols, var)
    v[inds]
  end

  # construct new table
  𝒯 = (; (var => genvar(var) for var in vars)...)
  newtab = 𝒯 |> Tables.materializer(tab)

  # new spatial data
  newgeotable = georef(newtab, tdom)

  newgeotable, nothing
end
