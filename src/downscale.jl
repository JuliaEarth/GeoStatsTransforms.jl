# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Downscale(f₁, f₂, ..., fₙ)

Downscale each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

Resulting values are obtained with the [`Transfer`](@ref) transform.

## Examples

```julia
Downscale(2, 2)
Downscale(3, 3, 2)
```
"""
struct Downscale{Dim} <: TableTransform
  factors::Dims{Dim}
end

Downscale(factors::Int...) = Downscale(factors)

function apply(transform::Downscale, geotable::AbstractGeoTable)
  gtb = geotable |> AbsoluteUnits()
  tab = values(gtb)
  grid = domain(gtb)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  # downscale the grid
  factors = _fitdims(transform.factors, paramdim(grid))
  tgrid = refine(grid, RegularRefinement(factors))

  # perform transfer
  pairs = map(vars) do var
    svals = Tables.getcolumn(cols, var)
    array = similar(svals, size(tgrid))
    titer = TileIterator(axes(array), factors)
    for (sind, tinds) in enumerate(titer)
      array[tinds...] .= svals[sind]
    end
    tvals = vec(array)
    var => tvals
  end

  # construct new table
  newtab = (; pairs...) |> Tables.materializer(tab)

  # new spatial data
  newgeotable = georef(newtab, tgrid)

  newgeotable, nothing
end
