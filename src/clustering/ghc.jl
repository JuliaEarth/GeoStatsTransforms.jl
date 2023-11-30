# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary functions and variables
uniform(h; λ) = (h ≤ λ)
triangular(h; λ) = (h ≤ λ) * (λ - h)
epanechnikov(h; λ) = (h ≤ λ) * (λ^2 - h^2)
const KERNFUN = Dict(:uniform => uniform, :triangular => triangular, :epanechnikov => epanechnikov)

"""
    GHC(k, λ; kern=:epanechnikov, link=:ward, as=:CLUSTER)

A transform for partitioning geospatial data into `k` clusters 
according to a range `λ` using Geostatistical Hierarchical
Clustering (GHC). The larger the range the more connected
are nearby samples.

## Parameters

* `k`    - Approximate number of clusters
* `λ`    - Approximate range of kernel function
* `kern` - Kernel function (`:uniform`, `:triangular` or `:epanechnikov`)
* `link` - Linkage function (`:single`, `:average`, `:complete`, `:ward` or `:ward_presquared`)
* `as`   - Cluster column name

## References

* Fouedjio, F. 2016. [A hierarchical clustering method for multivariate geostatistical data]
  (https://www.sciencedirect.com/science/article/abs/pii/S2211675316300367)

### Notes

- The range parameter controls the sparsity pattern of the pairwise
  distances, which can greatly affect the computational performance
  of the GHC algorithm. We recommend choosing a range that is small
  enough to connect nearby samples. For example, clustering data over
  a 100x100 Cartesian grid with unit spacing is possible with `λ=1.0`
  or `λ=2.0` but the problem starts to become computationally unfeasible
  around `λ=10.0` due to the density of points.
"""
struct GHC <: ClusteringTransform
  k::Int
  λ::Float64
  kern::Symbol
  link::Symbol
  as::Symbol
end

function GHC(k, λ; kern=:epanechnikov, link=:ward, as=:CLUSTER)
  # sanity checks
  @assert k > 0 "invalid number of clusters"
  @assert λ > 0 "invalid kernel range"
  @assert kern ∈ [:uniform, :triangular, :epanechnikov] "invalid kernel function"
  @assert link ∈ [:single, :average, :complete, :ward, :ward_presquared] "invalid linkage function"
  GHC(k, λ, kern, link, Symbol(as))
end

function apply(transform::GHC, geotable)
  # GHC parameters
  k = transform.k
  λ = transform.λ
  kern = transform.kern
  link = transform.link

  # all covariates must be continuous
  values(geotable) |> SciTypeAssertion{Continuous}()

  # dissimilarity matrix
  D = ghc_dissimilarity_matrix(geotable, kern, λ)

  # classical hierarchical clustering
  tree = hclust(D, linkage=link)

  # cut tree to produce clusters
  labels = cutree(tree, k=k)

  newtable = (; transform.as => categorical(labels))
  newgeotable = georef(newtable, domain(geotable))

  newgeotable, nothing
end

function ghc_dissimilarity_matrix(geotable, kern, λ)
  # retrieve domain/table
  𝒟 = domain(geotable)
  𝒯 = values(geotable)

  # kernel matrix
  K = ghc_kernel_matrix(kern, λ, 𝒟)

  # difference matrices
  Δ = ghc_diff_matrices(𝒯)

  # sum of cross-variograms
  Γ = ghc_variogram_sum(K, Δ)
end

function ghc_kernel_matrix(kern, λ, 𝒟)
  # kernel function
  fn = KERNFUN[kern]
  Kλ(h) = fn(h, λ=λ)

  # collect coordinates
  coords = [coordinates(centroid(𝒟, i)) for i in 1:nelements(𝒟)]

  # lag matrix
  H = pairwise(Euclidean(), coords)

  # kernel matrix
  K = Kλ.(H)

  # return sparse version
  sparse(K)
end

function ghc_diff_matrices(𝒯)
  # features must be standardized
  𝒮 = ghc_standardize(𝒯)

  # retrieve standardized features
  cols = Tables.columns(𝒮)
  vars = Tables.columnnames(cols)

  # distance matrices
  D = map(vars) do var
    z = Tables.getcolumn(cols, var)
    pairwise(Euclidean(), z)
  end

  # number of covariates
  p = length(vars)

  # one matrix per covariate pair
  Δ = Matrix{Matrix{Float64}}(undef, p, p)
  @inbounds for j in 1:p
    for i in (j + 1):p
      Δ[i, j] = D[i] .* D[j]
    end
    Δ[j, j] = D[j] .* D[j]
    for i in 1:(j - 1)
      Δ[i, j] = Δ[j, i] # leverage the symmetry
    end
  end

  Δ
end

function ghc_variogram_sum(K, Δ)
  n = size(K, 1)
  Γ = zeros(n, n)
  for Δₒ in Δ # for each covariate pair
    # update lower triangular matrix
    @inbounds for j in 1:n
      kj = K[:, j]
      for i in (j + 1):n
        ki = K[:, i]
        Kij = kron(ki, kj)
        I, W = findnz(Kij)
        num = sum(W .* Δₒ[I], init=0.0)
        den = sum(W, init=0.0)
        iszero(den) || (Γ[i, j] += (1 / 2) * (num / den))
      end
    end
  end

  # mirror upper triangular matrix
  @inbounds for j in 1:n
    Γ[j, j] = 0.0
    for i in 1:(j - 1)
      Γ[i, j] = Γ[j, i] # leverage the symmetry
    end
  end

  Γ
end

function ghc_standardize(𝒯)
  cols = Tables.columns(𝒯)
  vars = Tables.columnnames(cols)
  zstd = map(vars) do var
    z = Tables.getcolumn(cols, var)
    μ = mean(z)
    σ = std(z, mean=μ)
    iszero(σ) ? zero(μ) : (z .- μ) ./ σ
  end
  (; zip(vars, zstd)...) |> Tables.materializer(𝒯)
end
