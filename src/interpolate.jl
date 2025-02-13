# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Interpolate(domain; [options])
  
Interpolate geotable on given `domain` (or vector of geometries)
using a set of `options`.

## Options

* `model` - Model from GeoStatsModels.jl (default to `NN()`)
* `point` - Perform interpolation on point support (default to `true`)
* `prob`  - Perform probabilistic interpolation (default to `false`)

## Examples

```julia
# nearest neighbor model
geotable |> Interpolate(grid, model=NN())

# inverse distance weighting model
geotable |> Interpolate(pset, model=IDW())
```

### Notes

The interpolation is performed with all the samples (i.e. rows of
the geotable) at once, which can be prohibitive for some models
like Kriging. If that is the case, prefer [`InterpolateNeighbors`](@ref).

See also [`InterpolateNeighbors`](@ref).
"""
struct Interpolate{D<:Domain,GM<:GeoStatsModel} <: TableTransform
  domain::D
  model::GM
  point::Bool
  prob::Bool
end

Interpolate(domain::Domain; model=NN(), point=true, prob=false) = Interpolate(domain, model, point, prob)

Interpolate(geoms::AbstractVector{<:Geometry}; kwargs...) = Interpolate(GeometrySet(geoms); kwargs...)

isrevertible(::Type{<:Interpolate}) = false

function apply(transform::Interpolate, geotable::AbstractGeoTable)
  interp = fitpredict(
    # forward arguments
    transform.model,
    geotable |> AbsoluteUnits(), # handle affine units
    transform.domain;
    point=transform.point,
    prob=transform.prob,
    neighbors=false
  )

  interp, nothing
end
