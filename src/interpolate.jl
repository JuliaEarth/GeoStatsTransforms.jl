# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Interpolate(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
    Interpolate([g₁, g₂, ..., gₙ], vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` or set of geometries `g₁`, `g₂`, ..., `gₙ`,
using geostatistical models `model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.

    Interpolate(domain, model=IDW(); [parameters])
    Interpolate([g₁, g₂, ..., gₙ], model=IDW(); [parameters])
  
Interpolate geospatial data on given `domain` or set of geometries `g₁`, `g₂`, ..., `gₙ`,
using geostatistical `model` for all variables.

## Parameters

* `point` - Perform interpolation on point support (default to `true`)
* `prob`  - Perform probabilistic interpolation (default to `false`)

See also [`InterpolateNeighbors`](@ref).
"""
struct Interpolate{D<:Domain} <: TableTransform
  domain::D
  colspecs::Vector{ColSpec}
  models::Vector{GeoStatsModel}
  point::Bool
  prob::Bool
end

Interpolate(domain::Domain, colspecs, models; point=true, prob=false) =
  Interpolate(domain, collect(ColSpec, colspecs), collect(GeoStatsModel, models), point, prob)

Interpolate(geoms::AbstractVector{<:Geometry}, colspecs, models; kwargs...) =
  Interpolate(GeometrySet(geoms), colspecs, models; kwargs...)

Interpolate(domain, model::GeoStatsModel=IDW(); kwargs...) =
  Interpolate(domain, [AllSpec()], [model]; kwargs...)

Interpolate(domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  Interpolate(domain, colspec.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:Interpolate}) = false

function apply(transform::Interpolate, geotable::AbstractGeoTable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  idom = transform.domain
  colspecs = transform.colspecs
  models = transform.models
  point = transform.point
  prob = transform.prob

  data = if point
    pset = PointSet(centroid(dom, i) for i in 1:nelements(dom))
    georef(values(geotable), pset) |> uadjust
  else
    uadjust(geotable)
  end

  # preprocess variable models
  varmodels = mapreduce(vcat, colspecs, models) do colspec, model
    fmodel = fit(model, data)
    svars = choose(colspec, vars)
    [var => fmodel for var in svars]
  end

  # prediction order
  inds = traverse(idom, LinearPath())

  # predict variable values
  function pred(var, fmodel)
    map(inds) do ind
      geom = point ? centroid(idom, ind) : idom[ind]
      pfun = prob ? predictprob : predict
      pfun(fmodel, var, geom)
    end
  end

  pairs = (var => pred(var, fmodel) for (var, fmodel) in varmodels)
  newtab = (; pairs...) |> Tables.materializer(tab)

  newgeotable = georef(newtab, idom)

  newgeotable, nothing
end
