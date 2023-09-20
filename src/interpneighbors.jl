# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InterpolateNeighbors(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` using geostatistical models
`model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.

    InterpolateNeighbors(domain, model=IDW(); [parameters])
  
Interpolate geospatial data on given `domain` using geostatistical `model` for all variables.

Unlike [`Interpolate`](@ref), this transform uses neighbor search methods to
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

See also [`Interpolate`](@ref).
"""
struct InterpolateNeighbors{D<:Domain,N,M} <: TableTransform
  domain::D
  colspecs::Vector{ColSpec}
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
  colspecs,
  models;
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  point=true,
  prob=false
) = InterpolateNeighbors(
  domain,
  collect(ColSpec, colspecs),
  collect(GeoStatsModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  point,
  prob
)

InterpolateNeighbors(domain::Domain, model::GeoStatsModel=IDW(); kwargs...) =
  InterpolateNeighbors(domain, [AllSpec()], [model]; kwargs...)

InterpolateNeighbors(domain::Domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
  InterpolateNeighbors(domain, colspec.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:InterpolateNeighbors}) = false

function apply(transform::InterpolateNeighbors, geotable::AbstractGeoTable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  idom = transform.domain
  colspecs = transform.colspecs
  models = transform.models
  minneighbors = transform.minneighbors
  maxneighbors = transform.maxneighbors
  neighborhood = transform.neighborhood
  distance = transform.distance
  point = transform.point
  prob = transform.prob

  nobs = nelements(dom)
  if maxneighbors > nobs || maxneighbors < 1
    @warn "Invalid maximum number of neighbors. Adjusting to $nobs..."
    maxneighbors = nobs
  end

  if minneighbors > maxneighbors || minneighbors < 1
    @warn "Invalid minimum number of neighbors. Adjusting to 1..."
    minneighbors = 1
  end

  data = if point
    pset = PointSet(centroid(dom, i) for i in 1:nobs)
    georef(values(geotable), pset) |> uadjust
  else
    uadjust(geotable)
  end

  # preprocess variable models
  varmodels = mapreduce(vcat, colspecs, models) do colspec, model
    svars = choose(colspec, vars)
    [var => model for var in svars]
  end

  # determine bounded search method
  searcher = searcher_ui(dom, maxneighbors, distance, neighborhood)

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # prediction order
  inds = traverse(idom, LinearPath())

  # predict variable values
  function pred(var, model)
    map(inds) do ind
      # centroid of estimation
      center = centroid(idom, ind)

      # find neighbors with data
      nneigh = search!(neighbors, center, searcher)

      # predict if enough neighbors
      if nneigh ≥ minneighbors
        # final set of neighbors
        ninds = view(neighbors, 1:nneigh)

        # view neighborhood with data
        samples = view(data, ninds)

        # fit model to data
        fmodel = fit(model, samples)

        # save prediction
        geom = point ? center : idom[ind]
        pfun = prob ? predictprob : predict
        pfun(fmodel, var, geom)
      else # missing prediction
        missing
      end
    end
  end

  pairs = (var => pred(var, model) for (var, model) in varmodels)
  newtab = (; pairs...) |> Tables.materializer(tab)

  newgeotable = georef(newtab, idom)

  newgeotable, nothing
end
