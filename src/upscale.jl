# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Upscale(f1, f2, ...)

Upscale each dimension of the grid by given factors `f1, f2, ...`.

This transform is equivalent to skipping entries of the grid
as in the pseudo-code `grid[1:f1:end, 1:f2:end, ...]`.

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

function _targetgrid(grid::CartesianGrid{Dim}, factors::Dims{Dim}) where {Dim}
  dims = size(grid) .÷ factors
  CartesianGrid(minimum(grid), maximum(grid); dims)
end

function _targetgrid(grid::RectilinearGrid{Dim}, factors::Dims{Dim}) where {Dim}
  xyz = Meshes.xyz(grid)
  dims = size(grid) .+ .!isperiodic(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  RectilinearGrid(ntuple(i -> xyz[i][ranges[i]], Dim))
end

function _targetgrid(grid::StructuredGrid{Dim}, factors::Dims{Dim}) where {Dim}
  XYZ = Meshes.XYZ(grid)
  dims = size(grid) .+ .!isperiodic(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  StructuredGrid(ntuple(i -> XYZ[i][ranges...], Dim))
end

function apply(transform::Upscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = _targetgrid(grid, transform.factors)
  newgeotable = geotable |> Aggregate(tgrid)
  newgeotable, nothing
end
