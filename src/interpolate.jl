# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Interpolate(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
    Interpolate([g₁, g₂, ..., gₙ], vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` or vector of geometries `[g₁, g₂, ..., gₙ]`,
using geostatistical models `model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.

    Interpolate(domain, model=NN(); [parameters])
    Interpolate([g₁, g₂, ..., gₙ], model=NN(); [parameters])
  
Interpolate geospatial data on given `domain` or vector of geometries `[g₁, g₂, ..., gₙ]`,
using geostatistical `model` for all variables.

## Parameters

* `point` - Perform interpolation on point support (default to `true`)
* `prob`  - Perform probabilistic interpolation (default to `false`)

See also [`InterpolateNeighbors`](@ref), [`InterpolateMissing`](@ref), [`InterpolateNaN`](@ref).
"""
struct Interpolate{D<:Domain} <: TableTransform
  domain::D
  selectors::Vector{ColumnSelector}
  models::Vector{GeoStatsModel}
  point::Bool
  prob::Bool
end

Interpolate(domain::Domain, selectors, models; point=true, prob=false) =
  Interpolate(domain, collect(ColumnSelector, selectors), collect(GeoStatsModel, models), point, prob)

Interpolate(geoms::AbstractVector{<:Geometry}, selectors, models; kwargs...) =
  Interpolate(GeometrySet(geoms), selectors, models; kwargs...)

Interpolate(domain, model::GeoStatsModel=NN(); kwargs...) = Interpolate(domain, [AllSelector()], [model]; kwargs...)

Interpolate(domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  Interpolate(domain, selector.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:Interpolate}) = false

function apply(transform::Interpolate, geotable::AbstractGeoTable)
  geotable′ = geotable |> AbsoluteUnits()

  cols = Tables.columns(values(geotable′))
  vars = Tables.columnnames(cols)

  dom = transform.domain
  point = transform.point
  prob = transform.prob

  selectors = transform.selectors
  models = transform.models

  interps = map(selectors, models) do selector, model
    gtb = geotable′[:, selector(vars)]
    fitpredict(model, gtb, dom; point, prob, neighbors=false)
  end

  newgeotable = reduce(hcat, interps)

  newgeotable, nothing
end
