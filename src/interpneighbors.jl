# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateNeighbors(domain; [parameters])
  
Interpolate geospatial data with neighbor search on given `domain`
(or vector of geometries) using a set of optional `parameters`.

## Parameters

* `model`        - Model from GeoStatsModels.jl (default to `NN()`)
* `path`         - Path over the domain (default to `LinearPath()`)
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

See also [`Interpolate`](@ref).
"""
struct InterpolateNeighbors{D<:Domain,GM<:GeoStatsModel,P,N,M} <: TableTransform
  domain::D
  model::GM
  path::P
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
  path=LinearPath(),
  point=true,
  prob=false,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean()
) = InterpolateNeighbors(domain, model, path, point, prob, minneighbors, maxneighbors, neighborhood, distance)

InterpolateNeighbors(geoms::AbstractVector{<:Geometry}; kwargs...) = InterpolateNeighbors(GeometrySet(geoms); kwargs...)

isrevertible(::Type{<:InterpolateNeighbors}) = false

function apply(transform::InterpolateNeighbors, geotable::AbstractGeoTable)
  interp = fitpredict(
    # forward arguments
    transform.model,
    geotable |> AbsoluteUnits(), # handle affine units
    transform.domain;
    path=transform.path,
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
