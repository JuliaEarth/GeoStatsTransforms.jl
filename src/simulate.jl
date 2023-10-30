# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
    Simulate(domain, vars₁ => process₁, ..., varsₙ => processₙ)
    Simulate(domain, nreals, vars₁ => process₁, ..., varsₙ => processₙ)
    Simulate([g₁, g₂, ..., gₙ], vars₁ => process₁, ..., varsₙ => processₙ)
    Simulate([g₁, g₂, ..., gₙ], nreals, vars₁ => process₁, ..., varsₙ => processₙ)

TODO
"""
struct Simulate{D<:Domain,R<:AbstractRNG} <: TableTransform 
  domain::D
  nreals::Int
  selectors::Vector{ColumnSelector}
  processes::Vector{GeoStatsProcess}
  rng::R
end

Simulate(domain::Domain, nreals::Int, selectors, processes, rng) = 
  Simulate(domain, nreals, collect(ColumnSelector, selectors), collect(GeoStatsProcess, processes), rng)

Simulate(geoms::AbstractVector{<:Geometry}, nreals::Int, selectors, processes, rng) =
  Simulate(GeometrySet(geoms), nreals, selectors, processes, rng)

Simulate(domain, nreals::Int, pairs::Pair{<:Any,<:GeoStatsProcess}...; rng=Random.default_rng()) =
  Simulate(domain, nreals, selector.(first.(pairs)), last.(pairs), rng)

Simulate(domain, pairs::Pair{<:Any,<:GeoStatsProcess}...; rng=Random.default_rng()) =
  Simulate(domain, 1, selector.(first.(pairs)), last.(pairs), rng)

isrevertible(::Type{<:Simulate}) = false

function apply(transform::Simulate, geotable::AbstractGeoTable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  (; domain, nreals, selectors, processes, rng) = transform
  ensembles = map(selectors, processes) do selector, process
    svars = selector(vars)
    data = geotable[:, svars]
    svars => rand(rng, process, domain, data, nreals)
  end

  pad = ndigits(nreals)
  newgeotable = mapreduce(hcat, ensembles) do (vars, ensemble)
    mapreduce(hcat, 1:nreals) do i
      gtb = ensemble[i]
      pairs = (v => "$(v)_$(string(i; pad))" for v in vars)
      gtb |> Rename(pairs...)
    end
  end

  newgeotable, nothing
end
