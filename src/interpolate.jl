# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Interpolate(
      domain,
      vars₁ => model₁, ..., varsₙ => modelₙ;
      minneighbors=nothing,
      maxneighbors=nothing,
      neighborhood=nothing,
      distance=Euclidean(),
      path=LinearPath()
    )
  
TODO
"""
struct Interpolate{D<:Domain,N,M,P} <: TableTransform
  domain::D
  colspecs::Vector{ColSpec}
  models::Vector{GeoStatsModel}
  minneighbors::Union{Int,Nothing}
  maxneighbors::Union{Int,Nothing}
  neighborhood::N
  distance::M
  path::P
end

Interpolate(
  domain::Domain,
  colspecs,
  models;
  minneighbors=nothing,
  maxneighbors=nothing,
  neighborhood=nothing,
  distance=Euclidean(),
  path=LinearPath()
) = Interpolate(
  domain,
  collect(ColSpec, colspecs),
  collect(GeoStatsModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  path
)

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

  nelms = nrow(geotable)
  nmin = isnothing(minneighbors) ? 1 : minneighbors
  nmax = isnothing(maxneighbors) ? nelms : maxneighbors

  if nmax > nelms || nmax < 1
    @warn "Invalid maximum number of neighbors. Adjusting to $nelms..."
    nmax = nelms
  end

  if nmin > nmax || nmin < 1
    @warn "Invalid minimum number of neighbors. Adjusting to 1..."
    nmin = 1
  end

  searcher = searcher_ui(dom, nmax, distance, neighborhood)

  # preprocess variable models
  varmodels = mapreduce(vcat, colspecs, models) do colspec, model
    svars = choose(colspec, vars)
    svars .=> model
  end

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, nmax)

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
      if nneigh < nmin
        missing
      else
        # final set of neighbors
        ninds = view(neighbors, 1:nneigh)
  
        # view neighborhood with data
        samples = view(geotable, ninds)
  
        # fit model to data
        fmodel = fit(model, samples)
  
        # save prediction
        predict(fmodel, var, idom[ind])
      end
    end
  end

  pairs = (var => pred(var, model) for (var, model) in varmodels)
  newtab = (; pairs...) |> Tables.materializer(tab)

  # new spatial data
  newgeotable = georef(newtab, idom)

  newgeotable, nothing
end
