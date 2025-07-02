# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Quenching(transiogram; [options])
  
Simulated quenching (Carle et al. 1998) with given
theoretical `transiogram` from GeoStatsFunctions.jl.

## Options

* `skip`         - Indices to skip during simulation (default to `[]`)
* `tol`          - Tolerance on relative error (default to `1e-2`)
* `maxiter`      - Maximum number of iterations (default to `10`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `rng`          - Random number generator (default to `Random.default_rng()`)

## Examples

```julia
# simulated quenching with Gaussian transiogram
Quenching(GaussianTransiogram())

# do not change value at first geometry of the domain
Quenching(SphericalTransiogram(), skip=[1])
```

## References

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct Quenching{T<:Transiogram,RNG} <: TableTransform
  transiogram::T
  skip::Vector{Int}
  tol::Float64
  maxiter::Int
  maxneighbors::Int
  rng::RNG
end

Quenching(func; skip=Int[], tol=1e-2, maxiter=10, maxneighbors=26, rng=Random.default_rng()) =
  Quenching(func, skip, tol, maxiter, maxneighbors, rng)

isrevertible(::Type{<:Quenching}) = false

function apply(transform::Quenching, geotable::AbstractGeoTable)
  # theoretical transiogram
  τ = transform.transiogram

  # sanity checks
  cols = Tables.columns(values(geotable))
  vars = Tables.columnnames(cols)
  vals = Tables.getcolumn(cols, first(vars))
  levs = levels(vals)
  isok = length(vars) == 1 && elscitype(vals) <: Categorical
  @assert isok "Quenching only defined for a single categorical variable"
  isok = length(levs) == nvariates(τ)
  @assert isok "Invalid transiogram for number of categorical values"

  # variable and domain
  var = first(vars)
  dom = domain(geotable)

  # auxiliary variables
  d = embeddim(dom)
  𝓁 = range(τ)
  ℒ = typeof(𝓁)

  # reference point
  p = Point(ntuple(i -> ℒ(0), d))

  # basis vectors
  v = ntuple(d) do j
    Vec(ntuple(i -> ℒ(i == j), d))
  end

  # indices where data can be changed
  nelm = nelements(dom)
  inds = setdiff(1:nelm, transform.skip)

  # efficient lookup of neighbors
  nmax = transform.maxneighbors
  ball = MetricBall(range(τ))
  searcher = KBallSearch(dom, nmax, ball)
  neighbors = Vector{Int}(undef, nmax)

  # objective function
  function objective(gtb)
    linds = _levelindices(levs, gtb)
    map(v) do vⱼ
      t = DirectionalTransiogram(vⱼ, gtb, var; maxlag=𝓁)
      hs = t.abscissas
      ts = t.ordinates
      map(enumerate(hs)) do (i, h)
        tmat = getindex.(ts, i)
        τmat = τ(p, p + ustrip(h) * vⱼ)
        t̂mat = view(τmat, linds, linds)
        sum(abs2, tmat .- t̂mat)
      end |> sum
    end |> sum
  end

  # main loop
  gtb = deepcopy(geotable)
  for _ in 1:(transform.maxiter)
    keepgoing = false
    @inbounds for ind in shuffle(transform.rng, inds)
      # search neighbors
      n = search!(neighbors, centroid(dom, ind), searcher)
      ninds = view(neighbors, 1:n)

      # skip if not enough neighbors
      n > 1 || continue

      # compute mode of neighbor values
      cols = Tables.columns(values(gtb))
      vals = Tables.getcolumn(cols, var)
      mode = _mode(levs, view(vals, ninds))

      # geotable with neighbors
      ntb = view(gtb, ninds)

      # objective before change
      val = vals[ind]
      obj = objective(ntb)

      # objective after change
      vals[ind] = mode
      newobj = objective(ntb)

      # keep going if relative error is large
      if abs(newobj - obj) / obj ≥ transform.tol
        keepgoing = true
      end

      # undo change if necessary
      newobj < obj || (vals[ind] = val)
    end
    keepgoing || break
  end

  gtb, nothing
end

function _levelindices(levs, gtb)
  cols = Tables.columns(values(gtb))
  vars = Tables.columnnames(cols)
  vals = Tables.getcolumn(cols, first(vars))
  indexin(levels(vals), levs)
end
