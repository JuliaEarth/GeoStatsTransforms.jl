# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Downscale(f₁, f₂, ..., fₙ)

Downscale each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

Resulting values are obtained with the [`Transfer`](@ref) transform.

# Examples

```julia
Downscale(2, 2)
Downscale(3, 3, 2)
```
"""
struct Downscale{Dim} <: TableTransform
  factors::Dims{Dim}
end

Downscale(factors::Int...) = Downscale(factors)

isrevertible(::Type{<:Downscale}) = false

function apply(transform::Downscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = refine(grid, RegularRefinement(transform.factors))
  newgeotable = geotable |> Transfer(tgrid)
  newgeotable, nothing
end
