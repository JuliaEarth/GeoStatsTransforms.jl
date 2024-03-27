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

function _downscale(grid::CartesianGrid{Dim}, factors::Dims{Dim}) where {Dim}
  dims = size(grid) .* factors
  CartesianGrid(minimum(grid), maximum(grid); dims)
end

function _downscale(grid::RectilinearGrid{Dim}, factors::Dims{Dim}) where {Dim}
  xyz = Meshes.xyz(grid)
  xyz′ = ntuple(i -> _refine(xyz[i], factors[i]), Dim)
  RectilinearGrid(xyz′)
end

function _downscale(grid::StructuredGrid{Dim}, factors::Dims{Dim}) where {Dim}
  XYZ′ = _XYZ(grid, factors)
  StructuredGrid(XYZ′)
end

function _refine(x, f)
  x′ = mapreduce(vcat, 1:(length(x) - 1)) do i
    range(x[i], x[i + 1], f + 1)[begin:(end - 1)]
  end
  push!(x′, last(x))
  x′
end

function _XYZ(grid::StructuredGrid{2}, factors::Dims{2})
  fᵢ, fⱼ = factors
  sᵢ, sⱼ = size(grid)
  us = 0:(1 / fᵢ):1
  vs = 0:(1 / fⱼ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))

  mat(quad) = [coordinates(quad(u, v)) for u in us, v in vs]
  M = [mat(grid[i, j]) for i in 1:sᵢ, j in 1:sⱼ]

  C = mapreduce(catⱼ, 1:sⱼ) do j
    Mⱼ = mapreduce(catᵢ, 1:sᵢ) do i
      Mᵢⱼ = M[i, j]
      i == sᵢ ? Mᵢⱼ : Mᵢⱼ[begin:(end - 1), :]
    end
    j == sⱼ ? Mⱼ : Mⱼ[:, begin:(end - 1)]
  end

  X = getindex.(C, 1)
  Y = getindex.(C, 2)
  (X, Y)
end

function _XYZ(grid::StructuredGrid{3}, factors::Dims{3})
  fᵢ, fⱼ, fₖ = factors
  sᵢ, sⱼ, sₖ = size(grid)
  us = 0:(1 / fᵢ):1
  vs = 0:(1 / fⱼ):1
  ws = 0:(1 / fₖ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))
  catₖ(A...) = cat(A..., dims=Val(3))

  mat(hex) = [coordinates(hex(u, v, w)) for u in us, v in vs, w in ws]
  M = [mat(grid[i, j, k]) for i in 1:sᵢ, j in 1:sⱼ, k in 1:sₖ]

  C = mapreduce(catₖ, 1:sₖ) do k
    Mₖ = mapreduce(catⱼ, 1:sⱼ) do j
      Mⱼₖ = mapreduce(catᵢ, 1:sᵢ) do i
        Mᵢⱼₖ = M[i, j, k]
        i == sᵢ ? Mᵢⱼₖ : Mᵢⱼₖ[begin:(end - 1), :, :]
      end
      j == sⱼ ? Mⱼₖ : Mⱼₖ[:, begin:(end - 1), :]
    end
    k == sₖ ? Mₖ : Mₖ[:, :, begin:(end - 1)]
  end

  X = getindex.(C, 1)
  Y = getindex.(C, 2)
  Z = getindex.(C, 3)
  (X, Y, Z)
end

function apply(transform::Downscale, geotable::AbstractGeoTable)
  grid = domain(geotable)
  tgrid = _downscale(grid, transform.factors)
  newgeotable = geotable |> Transfer(tgrid)
  newgeotable, nothing
end
