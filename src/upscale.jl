# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Upscale(scales...)

TODO

# Examples

```julia
Upscale(2, 2)
Upscale(3, 3, 2)
```
"""
struct Upscale{Dim,S<:ColumnSelector}
  scales::Dims{Dim}
end

Upscale(scales::Int...) = Upscale(scales)

isrevertible(::Type{<:Upscale}) = false

function _targetgrid(grid::CartesianGrid{Dim}, scales::Dims{Dim}) where {Dim}
  dims = size(grid) .÷ scales
  CartesianGrid(minimum(grid), maximum(grid); dims)
end

function _targetgrid(grid::RectilinearGrid{Dim}, scales::Dims{Dim}) where {Dim}
  xyz = Meshes.xyz(grid)
  dims = size(grid) .+ 1
  ranges = ntuple(i -> 1:scales[i]:dims[i], Dim)
  RectilinearGrid(ntuple(i -> xyz[i][ranges[i]], Dim))
end

function _targetgrid(grid::StructuredGrid{Dim}, scales::Dims{Dim}) where {Dim}
  XYZ = Meshes.XYZ(grid)
  dims = size(grid) .+ 1
  ranges = ntuple(i -> 1:scales[i]:dims[i], Dim)
  StructuredGrid(ntuple(i -> XYZ[i][ranges...], Dim))
end

function apply(transform::Upscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = _targetgrid(grid, transform.scales)
  newgeotable = geotable |> Aggregate(tgrid)
  newgeotable, nothing
end
