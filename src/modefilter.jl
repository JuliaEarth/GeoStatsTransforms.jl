# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ModeFilter(func; [options])
  
Replace categorical values with mode of neighboring values.

## Options

* `skip`         - Indices to skip in domain (default to `[]`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `rng`          - Random number generator (default to `Random.default_rng()`)
"""
struct ModeFilter{RNG} <: TableTransform
  skip::Vector{Int}
  maxneighbors::Int
  rng::RNG
end

ModeFilter(; skip=Int[], maxneighbors=26, rng=Random.default_rng()) = ModeFilter(skip, maxneighbors, rng)

isrevertible(::Type{<:ModeFilter}) = false

function apply(transform::ModeFilter, geotable::AbstractGeoTable)
  # sanity checks
  cols = Tables.columns(values(geotable))
  vars = Tables.columnnames(cols)
  vals = Tables.getcolumn(cols, first(vars))
  levs = levels(vals)
  isok = length(vars) == 1 && elscitype(vals) <: Categorical
  @assert isok "ModeFilter only defined for a single categorical variable"

  # variable and domain
  var = first(vars)
  dom = domain(geotable)

  # efficient lookup of neighbors
  nelm = nelements(dom)
  nmax = transform.maxneighbors
  searcher = KNearestSearch(dom, nmax)
  neighbors = Vector{Int}(undef, nmax)

  # indices where data can be changed
  cinds = setdiff(1:nelm, transform.skip)

  # initialize buffer with results
  newvals = deepcopy(vals)

  # initialize mode filter mask
  mask = trues(nelm)

  # main mode filter loop
  @inbounds for ind in shuffle(transform.rng, cinds)
    # disable current index
    mask[ind] = false

    # center of target location
    center = centroid(dom, ind)

    # search neighbors with disabled index
    n = search!(neighbors, center, searcher, mask=mask)

    # skip if not enough neighbors
    n > 1 || continue

    # neighboring values
    inds = view(neighbors, 1:n)
    vals = view(newvals, inds)

    # update index with mode
    newvals[ind] = _mode(levs, vals)

    # re-enable current index and continue
    mask[ind] = true
  end

  # georeference results
  newgeotable = georef((; var => newvals), dom)

  newgeotable, nothing
end
