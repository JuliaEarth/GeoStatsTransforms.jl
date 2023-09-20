# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateGlobal(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` using geostatistical models
`model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.
Unlike [`Interpolate`](@ref), this transform performs the `fit` of the model once
with all the available data instead of multiple `fit` calls with neighborhoods.

## Parameters

* `path`  - The path algorithm used to iterate over the domain (default to `LinearPath()`)
* `point` - Perform interpolation on point support (default to `true`)
* `prob`  - Perform probabilistic interpolation (default to `false`)

See also [`Interpolate`](@ref).
"""
struct InterpolateGlobal{D<:Domain,P} <: TableTransform
  domain::D
  colspecs::Vector{ColSpec}
  models::Vector{GeoStatsModel}
  path::P
  point::Bool
  prob::Bool
end

InterpolateGlobal(domain::Domain, colspecs, models; path=LinearPath(), point=true, prob=false) =
  InterpolateGlobal(domain, collect(ColSpec, colspecs), collect(GeoStatsModel, models), path, point, prob)

InterpolateGlobal(domain::Domain; kwargs...) = InterpolateGlobal(domain, [AllSpec()], [IDW()]; kwargs...)

InterpolateGlobal(domain::Domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  InterpolateGlobal(domain, colspec.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:InterpolateGlobal}) = false

function apply(transform::InterpolateGlobal, geotable::AbstractGeoTable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  idom = transform.domain
  colspecs = transform.colspecs
  models = transform.models
  path = transform.path
  point = transform.point
  prob = transform.prob

  data = if point
    pset = PointSet(centroid(dom, i) for i in 1:nelements(dom))
    georef(values(geotable), pset)
  else
    geotable
  end

  # preprocess variable models
  varmodels = mapreduce(vcat, colspecs, models) do colspec, model
    fmodel = fit(model, data)
    svars = choose(colspec, vars)
    [var => fmodel for var in svars]
  end

  # prediction order
  inds = traverse(idom, path)

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
