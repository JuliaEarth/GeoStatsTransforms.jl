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

function _upscale(grid::CartesianGrid{Dim}, factors::Dims{Dim}) where {Dim}
  dims = size(grid) .÷ factors
  CartesianGrid(minimum(grid), maximum(grid); dims)
end

function _upscale(grid::RectilinearGrid{Dim}, factors::Dims{Dim}) where {Dim}
  xyz = Meshes.xyz(grid)
  dims = size(grid) .+ .!isperiodic(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  xyz′ = ntuple(i -> xyz[i][ranges[i]], Dim)
  RectilinearGrid(xyz′)
end

function _upscale(grid::StructuredGrid{Dim}, factors::Dims{Dim}) where {Dim}
  XYZ = Meshes.XYZ(grid)
  dims = size(grid) .+ .!isperiodic(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  XYZ′ = ntuple(i -> XYZ[i][ranges...], Dim)
  StructuredGrid(XYZ′)
end

function apply(transform::Upscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = _upscale(grid, transform.factors)
  newgeotable = geotable |> Aggregate(tgrid)
  newgeotable, nothing
end
