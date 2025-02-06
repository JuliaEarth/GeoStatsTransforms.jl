# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary functions and variables
uniform(h; λ) = (h ≤ λ)
triangular(h; λ) = (h ≤ λ) * (λ - h)
epanechnikov(h; λ) = (h ≤ λ) * (λ^2 - h^2)
const KERNFUN = Dict(:uniform => uniform, :triangular => triangular, :epanechnikov => epanechnikov)

"""
    GHC(k, λ; nmax=2000, kern=:epanechnikov, link=:ward, as=:cluster)

A transform for partitioning geospatial data into `k` clusters 
according to a range `λ` using Geostatistical Hierarchical
Clustering (GHC). The larger the range the more connected
are nearby samples.

## Parameters

* `k`    - Approximate number of clusters
* `λ`    - Approximate range of kernel function in length units
* `nmax` - Maximum number of observations to use in dissimilarity matrix
* `kern` - Kernel function (`:uniform`, `:triangular` or `:epanechnikov`)
* `link` - Linkage function (`:single`, `:average`, `:complete`, `:ward` or `:ward_presquared`)
* `as`   - Variable name used to store clustering results

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
struct GHC{ℒ<:Len} <: ClusteringTransform
  k::Int
  λ::ℒ
  nmax::Int
  kern::Symbol
  link::Symbol
  as::Symbol
  GHC(k, λ::ℒ, nmax, kern, link, as) where {ℒ<:Len} = new{float(ℒ)}(k, λ, nmax, kern, link, as)
end

function GHC(k, λ::Len; nmax=2000, kern=:epanechnikov, link=:ward, as=:cluster)
  # sanity checks
  @assert k > 0 "number of cluster must be positive"
  @assert λ > zero(λ) "kernel range must be positive"
  @assert nmax > 0 "maximum number of observations must be positive"
  @assert kern ∈ [:uniform, :triangular, :epanechnikov] "invalid kernel function"
  @assert link ∈ [:single, :average, :complete, :ward, :ward_presquared] "invalid linkage function"
  GHC(k, λ, nmax, kern, link, Symbol(as))
end

GHC(k, λ; kwargs...) = GHC(k, _addunit(λ, u"m"); kwargs...)

function apply(transform::GHC, geotable)
  # GHC parameters
  k = transform.k
  λ = transform.λ
  nmax = transform.nmax
  kern = transform.kern
  link = transform.link

  # all covariates must be continuous
  values(geotable) |> Assert(cond=x -> elscitype(x) <: Continuous)

  # sub-sample geotable to fit dissimilarity matrix
  gtb = ghc_subsample(geotable, nmax)

  # dissimilarity matrix
  D = ghc_dissimilarity_matrix(gtb, kern, λ)

  # classical hierarchical clustering
  tree = hclust(D, linkage=link)

  # cut tree to produce clusters
  labels = cutree(tree, k=k)

  # georeference categorical labels
  newtab = (; transform.as => categorical(labels))
  newgtb = georef(newtab, domain(gtb))

  # interpolate neighbors in case of sub-sampling
  newgeotable = if nrow(newgtb) < nrow(geotable)
    newgtb |> InterpolateNeighbors(domain(geotable), model=NN())
  else
    newgtb
  end

  newgeotable, nothing
end

function ghc_subsample(geotable, nmax)
  nobs = nrow(geotable)
  tran = nobs > nmax ? Sample(nmax, replace=false) : Identity()
  geotable |> tran
end

function ghc_dissimilarity_matrix(geotable, kern, λ)
  # retrieve domain/table
  𝒟 = domain(geotable)
  𝒯 = values(geotable)

  # kernel matrix
  K = ghc_kern_matrix(kern, λ, 𝒟)

  # features must be standardized
  𝒮 = ghc_standardize(𝒯)

  # retrieve feature columns
  cols = Tables.columns(𝒮)
  vars = Tables.columnnames(cols)

  # number of covariates
  p = length(vars)

  # number of observations
  n = size(K, 1)

  # dissimilarity matrix
  D = zeros(n, n)
  @inbounds for j in 1:p # for each pair of covariates
    zⱼ = Tables.getcolumn(cols, j)
    for i in j:p
      zᵢ = Tables.getcolumn(cols, i)

      # difference matrix for covariate pair
      Δ = ghc_diff_matrix(zᵢ, zⱼ)

      # contribution to dissimilarity matrix
      for l in 1:n
        Kₗ = K[:, l]
        for k in (l + 1):n
          Kₖ = K[:, k]
          Kₖₗ = kron(Kₗ, Kₖ) # faster Kₖ * transpose(Kₗ)
          I, W = findnz(Kₖₗ)
          num = sum(W .* Δ[I], init=zero(eltype(W)))
          den = sum(W, init=zero(eltype(W)))
          iszero(den) || (D[k, l] += (1 / 2) * (num / den))
        end
        D[l, l] = 0.0
        for k in 1:(l - 1)
          D[k, l] = D[l, k] # leverage symmetry
        end
      end
    end
  end

  D
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

function ghc_kern_matrix(kern, λ, 𝒟)
  # kernel function
  fn = KERNFUN[kern]
  Kλ(h) = fn(h, λ=λ)

  # collect coordinates
  coords = [to(centroid(𝒟, i)) for i in 1:nelements(𝒟)]

  # lag matrix
  H = pairwise(Euclidean(), coords)

  # kernel matrix
  K = ustrip.(Kλ.(H))

  # return sparse matrix
  sparse(K)
end

function ghc_diff_matrix(zᵢ, zⱼ)
  n = length(zᵢ)
  Δ = zeros(n, n)
  @inbounds for l in 1:n
    for k in (l + 1):n
      Δ[k, l] = (zᵢ[k] - zᵢ[l]) * (zⱼ[k] - zⱼ[l])
    end
    Δ[l, l] = 0.0
    for k in 1:(l - 1)
      Δ[k, l] = Δ[l, k] # leverage symmetry
    end
  end
  Δ
end
