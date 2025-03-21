# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Quenching(transiogram; [options])
  
Simulated quenching (Carle et al. 1998) with given
theoretical `transiogram` from GeoStatsFunctions.jl.

## Options

* `tol`     - Tolerance on relative error (default to `1e-2`)
* `maxiter` - Maximum number of iterations (default to `10`)
* `skip`    - Indices to skip during simulation (default to `[]`)

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
struct Quenching{T<:Transiogram} <: TableTransform
  transiogram::T
  tol::Float64
  maxiter::Int
  skip::Vector{Int}
end

Quenching(func; tol=1e-2, maxiter=10, skip=Int[]) = Quenching(func, tol, maxiter, skip)

isrevertible(::Type{<:Quenching}) = false

function apply(transform::Quenching, geotable::AbstractGeoTable)
  # theoretical transiogram
  τ = transform.transiogram

  # sanity checks
  cols = Tables.columns(values(geotable))
  vars = Tables.columnnames(cols)
  vals = Tables.getcolumn(cols, first(vars))
  isok = length(vars) == 1 && elscitype(vals) <: Categorical
  @assert isok "Quenching only defined for a single categorical variable"
  isok = length(levels(vals)) == nvariates(τ)
  @assert isok "Invalid transiogram for number of categorical values"

  # variable and domain
  var = first(vars)
  dom = domain(geotable)

  # indices where data can be changed
  nelm = nelements(dom)
  inds = setdiff(1:nelm, transform.skip)

  # searcher for efficient lookup of neighbors
  nmax = 26
  searcher = KNearestSearch(dom, nmax)
  neighbors = Vector{Int}(undef, nmax)

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

  # objective function
  function objective(gtb)
    map(v) do vⱼ
      t = DirectionalTransiogram(vⱼ, gtb, var; nlags=5, maxlag=𝓁)
      hs = t.abscissas
      ts = t.ordinates
      map(enumerate(hs)) do (i, h)
        tmat = getindex.(ts, i)
        τmat = τ(p, p + ustrip(h) * vⱼ)
        sum(abs2, tmat .- τmat)
      end |> sum
    end |> sum
  end

  # quenching step
  function quenching(gtb)
    cols = Tables.columns(values(gtb))
    vals = Tables.getcolumn(cols, var)
    @inbounds for i in shuffle(inds)
      c = centroid(dom, i)
      n = search!(neighbors, c, searcher)
      js = view(neighbors, 1:n)
      vs = view(vals, js)
      vals[i] = _mode(vs)
    end
    georef((; var => vals), dom)
  end

  # main loop
  iter = 0
  gtb = deepcopy(geotable)
  obj = objective(gtb)
  while iter < transform.maxiter
    # quenching step
    newgtb = quenching(gtb)
    newobj = objective(newgtb)

    # relative error
    error = abs((newobj - obj) / obj)

    # update in case of improvement
    if newobj < obj
      gtb = newgtb
      obj = newobj
    end

    # break in case of low error
    error < transform.tol && break

    iter += 1
  end

  gtb, nothing
end

function _mode(vs)
  c = Dict(unique(vs) .=> 0)
  for v in vs
    c[v] += 1
  end
  argmax(c)
end
