# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaxPosterior(func; [options])
  
Replace categorical values with maximum a posteriori
values given neighboring values and a theoretical
`func`tion from GeoStatsFunctions.jl.

## Options

* `skip`         - Indices to skip in domain (default to `[]`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `rng`          - Random number generator (default to `Random.default_rng()`)

## Examples

```julia
# maximum a posteriori with Gaussian transiogram
MaxPosterior(GaussianTransiogram())

# do not change value at first geometry of the domain
MaxPosterior(SphericalTransiogram(), skip=[1])
```

## References

* Deutsch, C. V. 2006. [A sequential indicator simulation
  program for categorical variables with point and block data:
  BlockSIS](https://www.sciencedirect.com/science/article/abs/pii/S0098300406000598)
"""
struct MaxPosterior{F<:GeoStatsFunction,RNG} <: TableTransform
  func::F
  skip::Vector{Int}
  maxneighbors::Int
  rng::RNG
end

MaxPosterior(func; skip=Int[], maxneighbors=26, rng=Random.default_rng()) = MaxPosterior(func, skip, maxneighbors, rng)

function apply(transform::MaxPosterior, geotable::AbstractGeoTable)
  # theoretical function
  f = transform.func

  # sanity checks
  cols = Tables.columns(values(geotable))
  vars = Tables.columnnames(cols)
  vals = Tables.getcolumn(cols, first(vars))
  levs = levels(vals)
  isok = length(vars) == 1 && elscitype(vals) <: Categorical
  @assert isok "MaxPosterior only defined for a single categorical variable"
  isok = length(levs) == nvariates(f)
  @assert isok "Invalid geostatistical function for number of categorical values"

  # variable and domain
  var = first(vars)
  dom = domain(geotable)

  # indicator variables in matrix form
  onehot = geotable |> OneHot(var) |> values
  onemat = stack(Tables.rowtable(onehot))

  # names of indicator variables
  ocols = Tables.columns(onehot)
  ovars = Tables.columnnames(ocols)

  # CoKriging model
  model = Kriging(f, _proportions(vals))

  # efficient lookup of neighbors
  nelm = nelements(dom)
  nmax = transform.maxneighbors
  ball = MetricBall(range(f))
  searcher = KBallSearch(dom, nmax, ball)
  neighbors = Vector{Int}(undef, nmax)

  # indices where data can be changed
  cinds = setdiff(1:nelm, transform.skip)

  # initialize buffer with results
  newvals = deepcopy(vals)

  # initialize maximum a posteriori mask
  mask = trues(nelm)

  # maximum a posteriori loop
  @inbounds for ind in shuffle(transform.rng, cinds)
    # disable current index
    mask[ind] = false

    # center of target location
    center = centroid(dom, ind)

    # search neighbors with disabled index
    n = search!(neighbors, center, searcher, mask=mask)

    # skip if not enough neighbors
    n > 1 || continue

    # neighborhood with data
    neigh = let
      inds = view(neighbors, 1:n)
      ndom = view(dom, inds)
      nmat = view(onemat, :, inds)
      ntup = ntuple(i -> view(nmat, i, :), length(ovars))
      ntab = NamedTuple{ovars}(ntup)
      georef(ntab, ndom)
    end

    # fit probability model
    fitted = GeoStatsModels.fit(model, neigh)

    # update index with posterior value
    if GeoStatsModels.status(fitted)
      prob = _predictprob(fitted, ovars, center)
      newvals[ind] = levs[argmax(prob)]
    end

    # re-enable current index and continue
    mask[ind] = true
  end

  # georeference results
  newgeotable = georef((; var => newvals), dom)

  newgeotable, nothing
end

function _proportions(vals)
  levs = levels(vals)
  counts = Dict{eltype(levs),Int}()
  for v in vals
    if v in keys(counts)
      counts[v] += 1
    else
      counts[v] = 1
    end
  end
  [counts[l] / length(vals) for l in levs]
end

function _predictprob(fitted, ovars, geom)
  p = GeoStatsModels.predict(fitted, ovars, geom)
  normalize(clamp.(p, 0, 1), 1)
end
