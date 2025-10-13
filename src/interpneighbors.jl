# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateNeighbors(domain; [options])
  
Interpolate geotable with neighbor search on given `domain`
(or vector of geometries) using a set of `options`.

## Options

* `model`        - Model from GeoStatsModels.jl (default to `NN()`)
* `point`        - Perform interpolation on point support (default to `true`)
* `prob`         - Perform probabilistic interpolation (default to `false`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance defined in Distances.jl (default to `Euclidean()`)

Two neighbor search methods are available:

* If a `neighborhood` is provided, local prediction is performed 
  by sliding the `neighborhood` in the domain (e.g. `MetricBall`).

* If a `neighborhood` is not provided, the prediction is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

## Examples

```julia
# polynomial model with 10 nearby samples
InterpolateNeighbors(grid, model=Polynomial(), maxneighbors=10)

# inverse distance weighting model with samples inside 100m radius
InterpolateNeighbors(pset, model=IDW(), neighborhood=MetricBall(100u"m"))
```

### Notes

The interpolation is performed with a subset of nearby samples.
This can lead to undesired artifacts, specially when the number
of samples is small. If that is that case, prefer [`Interpolate`](@ref).

See also [`Interpolate`](@ref).
"""
struct InterpolateNeighbors{D<:Domain,GM<:GeoStatsModel,N,M} <: TableTransform
  domain::D
  model::GM
  point::Bool
  prob::Bool
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
end

InterpolateNeighbors(
  domain::Domain;
  model=NN(),
  point=true,
  prob=false,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean()
) = InterpolateNeighbors(domain, model, point, prob, minneighbors, maxneighbors, neighborhood, distance)

InterpolateNeighbors(geoms::AbstractVector{<:Geometry}; kwargs...) = InterpolateNeighbors(GeometrySet(geoms); kwargs...)

function apply(transform::InterpolateNeighbors, geotable::AbstractGeoTable)
  interp = fitpredict(
    # forward arguments
    transform.model,
    geotable |> AbsoluteUnits(),
    transform.domain;
    point=transform.point,
    prob=transform.prob,
    neighbors=true,
    minneighbors=transform.minneighbors,
    maxneighbors=transform.maxneighbors,
    neighborhood=transform.neighborhood,
    distance=transform.distance
  )

  interp, nothing
end
