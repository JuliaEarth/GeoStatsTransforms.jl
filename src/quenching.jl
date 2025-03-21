# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Quenching(transiogram; maxiter=100, skip=[])
  
Simulated quenching (Carle et al. 1998) with given
theoretical `transiogram` from GeoStatsFunctions.jl.

Optionally, specify a maximum number of iterations
`maxiter`, and `skip` indices in the domain during
the simulation (i.e., optimization).

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
  maxiter::Int
  skip::Vector{Int}
end

Quenching(func; maxiter=100, skip=Int[]) = Quenching(func, maxiter, skip)

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
      j = rand((i - 1, i + 1))
      while j < 1 || j > nelm
        j = rand((i - 1, i + 1))
      end
      vals[i] = vals[j]
    end
    georef((; var => vals), domain(gtb))
  end

  # main loop
  iter = 0
  gtb = deepcopy(geotable)
  obj = objective(gtb)
  while iter < transform.maxiter
    # quenching step
    newgtb = quenching(gtb)
    newobj = objective(newgtb)

    # update upon improvement
    if newobj < obj
      gtb = newgtb
      obj = newobj
    end

    iter += 1
  end

  gtb, nothing
end

function _proportions(vals)
  levs = levels(vals)
  count = Dict(levs .=> 0)
  for v in vals
    count[v] += 1
  end
  [count[l] for l in levs] ./ length(levs)
end
