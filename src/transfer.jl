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
  table = values(gtb)
  cols = Tables.columns(table)
  vars = Tables.columnnames(cols)

  # source and target domains
  sdom = domain(gtb)
  tdom = transform.domain

  # perform transfer
  newcols = _transfer(sdom, tdom, cols, vars)
  newtable = (; newcols...) |> Tables.materializer(table)

  georef(newtable, tdom), nothing
end

function _transfer(sdom, tdom, cols, vars)
  if sdom isa Grid && tdom isa Grid && extrema(sdom) == extrema(tdom) && all(iszero, size(tdom) .% size(sdom))
    # we have two grids overlaid, and can rely on
    # tiled iteration for efficient transfer
    _gridtransfer(sdom, tdom, cols, vars)
  else
    # general case with knn search
    _knntransfer(sdom, tdom, cols, vars)
  end
end

function _gridtransfer(sdom, tdom, cols, vars)
  # determine tile size for tiled iteration
  tilesize = size(tdom) .÷ size(sdom)
  if any(<(1), tilesize)
    # fallback to general case with knn search
    _knntransfer(sdom, tdom, cols, vars)
  else
    # perform transfer with tiled iteration
    map(vars) do var
      svals = Tables.getcolumn(cols, var)
      array = similar(svals, size(tdom))
      titer = TileIterator(axes(array), tilesize)
      for (sind, tinds) in enumerate(titer)
        array[tinds...] .= svals[sind]
      end
      tvals = vec(array)
      var => tvals
    end
  end
end

function _knntransfer(sdom, tdom, cols, vars)
  # find nearest elements in source domain
  knn = KNearestSearch(sdom, 1)
  near = tmap(1:nelements(tdom)) do i
    first(search(centroid(tdom, i), knn))
  end

  # perform transfer
  map(vars) do var
    svals = Tables.getcolumn(cols, var)
    tvals = svals[near]
    var => tvals
  end
end
