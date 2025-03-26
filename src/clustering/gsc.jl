# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GSC(k, m; σ=1.0, tol=1e-4, maxiter=10, weights=nothing, as=:cluster)

A transform for partitioning geospatial data into `k` clusters
using Geostatistical Spectral Clustering (GSC).

## Parameters

* `k`       - Desired number of clusters
* `m`       - Multiplicative factor for adjacent weights

## Options

* `σ`       - Standard deviation for exponential model (default to `1.0`)
* `tol`     - Tolerance of k-means algorithm (default to `1e-4`)
* `maxiter` - Maximum number of iterations (default to `10`)
* `weights` - Dictionary with weights for each attribute (default to `nothing`)
* `as`      - Variable name used to store clustering results

## References

* Romary et al. 2015. [Unsupervised classification of multivariate
  geostatistical data: Two algorithms]
  (https://www.sciencedirect.com/science/article/pii/S0098300415001314)

## Notes

- The algorithm implemented here is slightly different than the algorithm
described in Romary et al. 2015. Instead of setting Wᵢⱼ = 0 when i <-/-> j,
we simply magnify the weight by a multiplicative factor Wᵢⱼ *= m when i <--> j.
This leads to dense matrices but also better results in practice.
"""
struct GSC{W} <: TableTransform
  k::Int
  m::Float64
  σ::Float64
  tol::Float64
  maxiter::Int
  weights::W
  as::Symbol
end

function GSC(k, m; σ=1.0, tol=1e-4, maxiter=10, weights=nothing, as=:cluster)
  # sanity checks
  @assert k > 0 "number of clusters must be positive"
  @assert m > 0 "multiplicative factor must be positive"
  @assert σ > 0 "standard deviation must be positive"
  GSC(k, m, σ, tol, maxiter, weights, Symbol(as))
end

function apply(transform::GSC, geotable::AbstractGeoTable)
  # retrieve table and domain
  𝒯 = values(geotable)
  𝒟 = domain(geotable)

  # retrieve parameters
  k = transform.k
  m = transform.m
  σ = transform.σ
  tol = transform.tol
  maxiter = transform.maxiter
  weights = transform.weights

  # table distance
  td = TableDistance(normalize=false, weights=weights)

  # adjacency matrix
  A = adjacencymatrix(𝒟)

  # weight matrix
  Δ = pairwise(td, 𝒯)
  E = @. exp(-Δ / σ^2)
  E[findall(!iszero, A)] .*= m
  W = sparse(E)

  # degree matrix
  d = vec(sum(W, dims=2))
  D = Diagonal(d)

  # Laplace matrix
  L = D^(-1 / 2) * W * D^(-1 / 2)

  # solve eigenproblem
  S, _ = partialschur(L, nev=k)
  _, V = partialeigen(S)

  # k-means with eigenvectors
  result = kmeans(V', k, tol=tol, maxiter=maxiter)
  labels = assignments(result)

  newtable = (; transform.as => categorical(labels))
  newgeotable = georef(newtable, domain(geotable))

  newgeotable, nothing
end
