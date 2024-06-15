# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateMissing(vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
    InterpolateMissing(vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on its own domain, using geostatistical models `model₁`, ..., `modelₙ` 
and non-missing values of the variables `vars₁`, ..., `varsₙ`.

    InterpolateMissing(model=NN(); [parameters])
    InterpolateMissing(model=NN(); [parameters])
  
Interpolate geospatial data on its own domain, using geostatistical `model` and non-missing values of all variables.

Just like [`InterpolateNeighbors`](@ref), this transform uses neighbor search methods to
fit geostatistical models at each interpolation location with a reduced number
of measurements.

## Parameters

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance defined in Distances.jl (default to `Euclidean()`)
* `point`        - Perform interpolation on point support (default to `true`)
* `prob`         - Perform probabilistic interpolation (default to `false`)

The `maxneighbors` parameter can be used to perform interpolation with
a subset of measurements per prediction location. If `maxneighbors`
is not provided, then all measurements are used.

Two `neighborhood` search methods are available:

* If a `neighborhood` is provided, local prediction is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is not provided, the prediction is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

See also [`InterpolateNaN`](@ref), [`InterpolateNeighbors`](@ref), [`Interpolate`](@ref).
"""
struct InterpolateMissing{N,M} <: TableTransform
  selectors::Vector{ColumnSelector}
  models::Vector{GeoStatsModel}
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
  point::Bool
  prob::Bool
end

InterpolateMissing(
  selectors,
  models;
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  point=true,
  prob=false
) = InterpolateMissing(
  collect(ColumnSelector, selectors),
  collect(GeoStatsModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  point,
  prob
)

InterpolateMissing(model::GeoStatsModel=NN(); kwargs...) = InterpolateMissing([AllSelector()], [model]; kwargs...)

InterpolateMissing(pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  InterpolateMissing(selector.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:InterpolateMissing}) = false

function _interp(geotable, selectors, models, droptrans; kwargs...)
  tab = values(geotable)
  dom = domain(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  interps = map(selectors, models) do selector, model
    svars = selector(vars)
    mapreduce(hcat, svars) do var
      data = geotable[:, [var]] |> droptrans
      fitpredict(model, data, dom; kwargs...)
    end
  end

  reduce(hcat, interps)
end

function apply(transform::InterpolateMissing, geotable::AbstractGeoTable)
  selectors = transform.selectors
  models = transform.models
  kwargs = (
    minneighbors=transform.minneighbors,
    maxneighbors=transform.maxneighbors,
    neighborhood=transform.neighborhood,
    distance=transform.distance,
    point=transform.point,
    prob=transform.prob
  )

  newgeotable = _interp(geotable, selectors, models, DropMissing(); kwargs...)

  newgeotable, nothing
end
