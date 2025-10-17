# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SLIC(k, m; tol=1e-4, maxiter=10, weights=nothing)

Assign labels to rows of geotable using the Simple Linear
Iterative Clustering (SLIC) algorithm.

The algorithm produces approximately `k` clusters by combining
a geospatial distance `dₛ` and a distance between variables `dᵥ`.
The tradeoff is controlled with a hyperparameter `m` in an
additive model `dₜ = √(dᵥ² + m²(dₛ/s)²)` where `s` is the
average spacing between cluster centroids.

## Options

* `tol`     - Tolerance of k-means algorithm (default to `1e-4`)
* `maxiter` - Maximum number of iterations (default to `10`)
* `weights` - Dictionary with weights for each variable (default to `nothing`)

## Examples

```julia
SLIC(10, 0.01) # default options
SLIC(10, 0.01, tol=1e-2) # set tolerance
```

## References

* Achanta et al. 2011. [SLIC superpixels compared to state-of-the-art
  superpixel methods](https://ieeexplore.ieee.org/document/6205760)

### Notes

The variables (or features) are standardized with the `StdFeats`
transform prior to the core algorithm in order to facilitate the
choice of the parameter `m`.
"""
struct SLIC{W} <: TableTransform
  k::Int
  m::Float64
  tol::Float64
  maxiter::Int
  weights::W
end

function SLIC(k::Int, m::Real; tol=1e-4, maxiter=10, weights=nothing)
  @assert tol > 0 "tolerance must be positive"
  @assert maxiter > 0 "maximum number of iterations must be positive"
  SLIC{typeof(weights)}(k, m, tol, maxiter, weights)
end

function apply(transform::SLIC, geotable::AbstractGeoTable)
  # retrieve parameters
  k = transform.k
  m = transform.m
  tol = transform.tol
  maxiter = transform.maxiter
  weights = transform.weights

  # normalize attributes
  Ω = geotable |> StdFeats()
  𝒟 = domain(Ω)

  # initial spacing of clusters
  s = slic_spacing(𝒟, k)

  # initialize cluster centers
  centers = slic_initialization(𝒟, s)

  # initialize search method
  searcher = BallSearch(𝒟, MetricBall(maximum(s)))

  # define table distance
  td = TableDistance(normalize=false, weights=weights)

  # pre-allocate memory for label and distance
  labels = fill(0, nelements(𝒟))
  dists = fill(Inf, nelements(𝒟))

  # Lloyd's (a.k.a. k-means) algorithm
  iter = 0
  δcur = mean(dists)
  while iter < maxiter
    slic_assignment!(Ω, searcher, td, m, s, centers, labels, dists)
    slic_update!(Ω, centers, labels)

    # average distance to centers
    δnew = mean(dists)

    # break upon convergence
    abs(δnew - δcur) / δcur < tol && break

    # update and continue
    δcur = δnew
    iter += 1
  end

  orphans = findall(iszero, labels)
  if length(orphans) > 0
    assigned = findall(!iszero, labels)
    𝒟₀ = view(𝒟, assigned)
    csearcher = KNearestSearch(𝒟₀, 1)

    for orphan in orphans
      p = centroid(𝒟, orphan)
      i = search(p, csearcher)[1]
      labels[orphan] = labels[assigned[i]]
    end
  end

  newtable = (; label=labels)
  newgeotable = georef(newtable, domain(geotable))

  newgeotable, nothing
end

slic_spacing(𝒟, k) = slic_srecursion(k, sides(boundingbox(𝒟)))

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

function slic_assignment!(geotable, searcher, td, m, s, centers, labels, dists)
  sₘ = maximum(s)
  𝒟 = domain(geotable)
  for (k, cₖ) in enumerate(centers)
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
    dᵥ = pairwise(td, V, vₖ)

    # total distance
    dₜ = @. √(dᵥ^2 + m^2 * (dₛ / sₘ)^2)

    @inbounds for (i, ind) in enumerate(inds)
      if dₜ[i] < dists[ind]
        dists[ind] = dₜ[i]
        labels[ind] = k
      end
    end
  end
end

function slic_update!(geotable, centers, labels)
  𝒟 = domain(geotable)
  for k in eachindex(centers)
    inds = findall(isequal(k), labels)
    X = (to(centroid(𝒟, i)) for i in inds)
    xₖ = [mean(X)]
    dₛ = pairwise(Euclidean(), X, xₖ)
    @inbounds centers[k] = inds[argmin(vec(dₛ))]
  end
end
