# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SLIC(k, m; tol=1e-4, maxiter=10, weights=nothing, as=:cluster)

A transform for clustering geospatial data into approximately `k`
clusters using Simple Linear Iterative Clustering (SLIC).

The transform produces clusters of samples that are spatially
connected based on a distance `dₛ` and that, at the same time,
are similar in terms of `vars` with distance `dᵥ`. The tradeoff
is controlled with a hyperparameter `m` in an additive model
`dₜ = √(dᵥ² + m²(dₛ/s)²)`.

## Parameters

* `k`       - Approximate number of clusters
* `m`       - Hyperparameter of SLIC model

## Options

* `tol`     - Tolerance of k-means algorithm (default to `1e-4`)
* `maxiter` - Maximum number of iterations (default to `10`)
* `weights` - Dictionary with weights for each attribute (default to `nothing`)
* `as`      - Variable name used to store clustering results

## References

* Achanta et al. 2011. [SLIC superpixels compared to state-of-the-art
  superpixel methods](https://ieeexplore.ieee.org/document/6205760)
"""
struct SLIC{W} <: TableTransform
  k::Int
  m::Float64
  tol::Float64
  maxiter::Int
  weights::W
  as::Symbol
end

function SLIC(k::Int, m::Real; tol=1e-4, maxiter=10, weights=nothing, as=:cluster)
  @assert tol > 0 "tolerance must be positive"
  @assert maxiter > 0 "maximum number of iterations must be positive"
  SLIC{typeof(weights)}(k, m, tol, maxiter, weights, Symbol(as))
end

function apply(transform::SLIC, geotable::AbstractGeoTable)
  # retrieve parameters
  w = transform.weights
  m = transform.m

  # normalize attributes
  𝒯 = TableDistances.normalize(values(geotable))
  Ω = georef(first(𝒯), domain(geotable))
  𝒟 = domain(Ω)

  # initial spacing of clusters
  s = slic_spacing(𝒟, transform)

  # initialize cluster centers
  c = slic_initialization(𝒟, s)

  # ball neighborhood search
  searcher = BallSearch(𝒟, MetricBall(maximum(s)))

  # pre-allocate memory for label and distance
  l = fill(0, nelements(𝒟))
  d = fill(Inf, nelements(𝒟))

  # performance parameters
  tol = transform.tol
  maxiter = transform.maxiter

  # Lloyd's (a.k.a. k-means) algorithm
  err, iter = Inf, 0
  while err > tol && iter < maxiter
    o = copy(c)

    slic_assignment!(Ω, searcher, w, m, s, c, l, d)
    slic_update!(Ω, c, l)

    err = norm(c - o) / norm(o)
    iter += 1
  end

  orphans = findall(iszero, l)
  if length(orphans) > 0
    assigned = findall(!iszero, l)
    𝒟₀ = view(𝒟, assigned)
    csearcher = KNearestSearch(𝒟₀, 1)

    for orphan in orphans
      p = centroid(𝒟, orphan)
      i = search(p, csearcher)[1]
      l[orphan] = l[assigned[i]]
    end
  end

  newtable = (; transform.as => categorical(l))
  newgeotable = georef(newtable, domain(geotable))

  newgeotable, nothing
end

slic_spacing(𝒟, transform) = slic_srecursion(transform.k, sides(boundingbox(𝒟)))

# given the desired number of clusters and the sides of the bounding box
# of the domain, returns the spacing for each dimension recursively
function slic_srecursion(k, l)
  d = length(l)

  # base case
  d == 1 && return [l[1] / k]

  # compute the spacing for the j-th dimension
  j = argmax(l)
  kⱼ = ceil(Int, k^(1 / d))
  sⱼ = l[j] / kⱼ

  # update the new k and l
  kₙ = ceil(Int, k / kⱼ)
  lₙ = l[[1:(j - 1); (j + 1):d]]

  # then recursively compute the spacing for the remaining dimensions
  s = slic_srecursion(kₙ, lₙ)

  [s[begin:(j - 1)]; [sⱼ]; s[j:end]]
end

function slic_initialization(𝒟, s)
  # efficient neighbor search
  searcher = KNearestSearch(𝒟, 1)

  # bounding box properties
  bbox = boundingbox(𝒟)
  lo, up = to.(extrema(bbox))

  # cluster centers
  clusters = Vector{Int}()
  neighbor = Vector{Int}(undef, 1)
  ranges = [(l + sᵢ / 2):sᵢ:u for (l, sᵢ, u) in zip(lo, s, up)]
  for x in Iterators.product(ranges...)
    search!(neighbor, Point(x), searcher)
    push!(clusters, neighbor[1])
  end

  unique(clusters)
end

function slic_assignment!(geotable, searcher, w, m, s, c, l, d)
  sₘ = maximum(s)
  𝒟 = domain(geotable)
  for (k, cₖ) in enumerate(c)
    inds = search(centroid(𝒟, cₖ), searcher)

    # distance between coordinates
    X = (to(centroid(𝒟, i)) for i in inds)
    xₖ = [to(centroid(𝒟, cₖ))]
    dₛ = pairwise(Euclidean(), X, xₖ)

    # distance between variables
    𝒮ᵢ = view(geotable, inds)
    𝒮ₖ = view(geotable, [cₖ])
    V = values(𝒮ᵢ)
    vₖ = values(𝒮ₖ)
    dᵥ = pairwise(TableDistance(normalize=false, weights=w), V, vₖ)

    # total distance
    dₜ = @. √(dᵥ^2 + m^2 * (dₛ / sₘ)^2)

    @inbounds for (i, ind) in enumerate(inds)
      if dₜ[i] < d[ind]
        d[ind] = dₜ[i]
        l[ind] = k
      end
    end
  end
end

function slic_update!(geotable, c, l)
  𝒟 = domain(geotable)
  for k in eachindex(c)
    inds = findall(isequal(k), l)
    X = (to(centroid(𝒟, i)) for i in inds)
    xₖ = [mean(X)]
    dₛ = pairwise(Euclidean(), X, xₖ)
    @inbounds c[k] = inds[argmin(vec(dₛ))]
  end
end
