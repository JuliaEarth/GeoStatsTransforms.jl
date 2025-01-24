# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateNeighbors(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
    InterpolateNeighbors([g₁, g₂, ..., gₙ], vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` or set of geometries `g₁`, `g₂`, ..., `gₙ`,
using geostatistical models `model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.

    InterpolateNeighbors(domain, model=NN(); [parameters])
    InterpolateNeighbors([g₁, g₂, ..., gₙ], model=NN(); [parameters])
  
Interpolate geospatial data on given `domain` or set of geometries `g₁`, `g₂`, ..., `gₙ`,
using geostatistical `model` for all variables.

Unlike [`Interpolate`](@ref), this transform uses neighbor search methods to
fit geostatistical models at each interpolation location with a reduced number
of measurements.

## Parameters

* `point`        - Perform interpolation on point support (default to `true`)
* `prob`         - Perform probabilistic interpolation (default to `false`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance defined in Distances.jl (default to `Euclidean()`)

The `maxneighbors` parameter can be used to perform interpolation with
a subset of measurements per prediction location. If `maxneighbors`
is not provided, then all measurements are used.

Two `neighborhood` search methods are available:

* If a `neighborhood` is provided, local prediction is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is not provided, the prediction is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

See also [`Interpolate`](@ref), [`InterpolateMissing`](@ref), [`InterpolateNaN`](@ref).
"""
struct InterpolateNeighbors{D<:Domain,N,M} <: TableTransform
  domain::D
  selectors::Vector{ColumnSelector}
  models::Vector{GeoStatsModel}
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
  point::Bool
  prob::Bool
end

InterpolateNeighbors(
  domain::Domain,
  selectors,
  models;
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  point=true,
  prob=false
) = InterpolateNeighbors(
  domain,
  collect(ColumnSelector, selectors),
  collect(GeoStatsModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  point,
  prob
)

InterpolateNeighbors(geoms::AbstractVector{<:Geometry}, selectors, models; kwargs...) =
  InterpolateNeighbors(GeometrySet(geoms), selectors, models; kwargs...)

InterpolateNeighbors(domain, model::GeoStatsModel=NN(); kwargs...) =
  InterpolateNeighbors(domain, [AllSelector()], [model]; kwargs...)

InterpolateNeighbors(domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  InterpolateNeighbors(domain, selector.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:InterpolateNeighbors}) = false

function apply(transform::InterpolateNeighbors, geotable::AbstractGeoTable)
  geotable′ = geotable |> AbsoluteUnits()

  cols = Tables.columns(values(geotable′))
  vars = Tables.columnnames(cols)

  dom = transform.domain
  point = transform.point
  prob = transform.prob
  minneighbors = transform.minneighbors
  maxneighbors = transform.maxneighbors
  neighborhood = transform.neighborhood
  distance = transform.distance

  selectors = transform.selectors
  models = transform.models

  interps = map(selectors, models) do selector, model
    gtb = geotable′[:, selector(vars)]
    fitpredict(model, gtb, dom; point, prob, minneighbors, maxneighbors, neighborhood, distance)
  end

  newgeotable = reduce(hcat, interps)

  newgeotable, nothing
end
