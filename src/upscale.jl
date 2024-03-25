# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Upscale(factors...)

TODO

# Examples

```julia
Upscale(2, 2)
Upscale(3, 3, 2)
```
"""
struct Upscale{Dim,S<:ColumnSelector}
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
  dims = size(grid) .+ open(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  RectilinearGrid(ntuple(i -> xyz[i][ranges[i]], Dim))
end

function _targetgrid(grid::StructuredGrid{Dim}, factors::Dims{Dim}) where {Dim}
  XYZ = Meshes.XYZ(grid)
  dims = size(grid) .+ open(grid)
  ranges = ntuple(i -> 1:factors[i]:dims[i], Dim)
  StructuredGrid(ntuple(i -> XYZ[i][ranges...], Dim))
end

function apply(transform::Upscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = _targetgrid(grid, transform.factors)
  newgeotable = geotable |> Aggregate(tgrid)
  newgeotable, nothing
end
