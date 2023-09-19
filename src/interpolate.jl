# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Interpolate(domain, vars₁ => model₁, ..., varsₙ => modelₙ; [parameters])
  
Interpolate geospatial data on given `domain` using geostatistical models
`model₁`, ..., `modelₙ` for variables `vars₁`, ..., `varsₙ`.

## Parameters

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance defined in Distances.jl (default to `Euclidean()`)
* `path`         - The path algorithm used to iterate over the domain (default to `LinearPath()`)
* `point`        - Perform interpolation on point support (default to `true`)
* `prob`         - Perform probabilistic interpolation (default to `false`)

The `maxneighbors` option can be used to perform interpolation
with a subset of measurements per prediction location. 
If `maxneighbors` is not provided, then all measurements are used.

Two `neighborhood` search methods are available:

* If a `neighborhood` is provided, local prediction is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is not provided, the prediction is performed 
  using `maxneighbors` nearest neighbors according to `distance`.
"""
struct Interpolate{D<:Domain,N,M,P} <: TableTransform
  domain::D
  colspecs::Vector{ColSpec}
  models::Vector{GeoStatsModel}
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
  path::P
  point::Bool
  prob::Bool
end

Interpolate(
  domain::Domain,
  colspecs,
  models;
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  path=LinearPath(),
  point=true,
  prob=false
) = Interpolate(
  domain,
  collect(ColSpec, colspecs),
  collect(GeoStatsModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  path,
  point,
  prob
)

Interpolate(domain::Domain; distance=Euclidean(), kwargs...) =
  Interpolate(domain, [AllSpec()], [IDW(1, distance)]; distance, kwargs...)

Interpolate(domain::Domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
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
  minneighbors = transform.minneighbors
  maxneighbors = transform.maxneighbors
  neighborhood = transform.neighborhood
  distance = transform.distance
  path = transform.path
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

  pset = PointSet(centroid(dom, i) for i in 1:nobs)
  data = point ? georef(values(geotable), pset) : geotable

  # preprocess variable models
  varmodels = mapreduce(vcat, colspecs, models) do colspec, model
    svars = choose(colspec, vars)
    [var => model for var in svars]
  end

  # determine bounded search method
  searcher = searcher_ui(pset, maxneighbors, distance, neighborhood)

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # prediction order
  inds = traverse(idom, path)

  # predict variable values
  function pred(var, model)
    map(inds) do ind
      # centroid of estimation
      center = centroid(idom, ind)

      # find neighbors with data
      nneigh = search!(neighbors, center, searcher)

      # skip if there are too few neighbors
      if nneigh < minneighbors
        missing
      else
        # final set of neighbors
        ninds = view(neighbors, 1:nneigh)

        # view neighborhood with data
        samples = view(data, ninds)

        # fit model to data
        fmodel = fit(model, samples)

        # save prediction
        geom = point ? center : dom[ind]
        pfun = prob ? predictprob : predict
        pfun(fmodel, var, geom)
      end
    end
  end

  pairs = (var => pred(var, model) for (var, model) in varmodels)
  newtab = (; pairs...) |> Tables.materializer(tab)

  newgeotable = georef(newtab, idom)

  newgeotable, nothing
end
