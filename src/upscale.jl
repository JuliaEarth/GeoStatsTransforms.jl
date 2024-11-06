# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Upscale(f₁, f₂, ..., fₙ)

Upscale each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

This transform is equivalent to skipping entries of the grid
as in the pseudo-code `grid[1:f₁:end, 1:f₂:end, ..., 1:fₙ:end]`.

Resulting values are obtained with the [`Aggregate`](@ref) transform
and its default aggregation functions.

# Examples

```julia
Upscale(2, 2)
Upscale(3, 3, 2)
```
"""
struct Upscale{Dim} <: TableTransform
  factors::Dims{Dim}
end

Upscale(factors::Int...) = Upscale(factors)

isrevertible(::Type{<:Upscale}) = false

function apply(transform::Upscale, geotable::AbstractGeoTable)
  gtb = _adjustunits(geotable)
  tab = values(gtb)
  grid = domain(gtb)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  # upscale the grid
  factors = _fitdims(transform.factors, paramdim(grid))
  tgrid = coarsen(grid, RegularCoarsening(factors))

  # aggregate columns
  pairs = map(vars) do var
    svals = Tables.getcolumn(cols, var)
    aggfun = _defaultagg(svals)
    array = reshape(svals, size(grid))
    titer = TileIterator(axes(array), factors)
    tvals = tmap(titer) do sinds
      aggfun(array[sinds...])
    end |> vec
    var => tvals
  end

  # construct new table
  newtab = (; pairs...) |> Tables.materializer(tab)

  # new spatial data
  newgeotable = georef(newtab, tgrid)

  newgeotable, nothing
end
