# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
    Simulate(domain, vars₁ => process₁, ..., varsₙ => processₙ; [parameters])
    Simulate(domain, nreals, vars₁ => process₁, ..., varsₙ => processₙ; [parameters])
    Simulate([g₁, g₂, ..., gₙ], vars₁ => process₁, ..., varsₙ => processₙ; [parameters])
    Simulate([g₁, g₂, ..., gₙ], nreals, vars₁ => process₁, ..., varsₙ => processₙ; [parameters])

Simulate `nreals` realizations of variables `varsᵢ` with geostatistical process
`processᵢ` over given `domain` or vector of geometries `[g₁, g₂, ..., gₙ]`.

## Parameters

* `rng`      - Random number generator (default to `Random.default_rng()`)
* `pool`     - Pool of worker processes (default to `[myid()]`)
* `threads`  - Number of threads (default to `cpucores()`)
* `progress` - Defines whether progress bar will be displayed (default to `true`)
"""
struct Simulate{D<:Domain,R<:AbstractRNG,K} <: TableTransform
  domain::D
  nreals::Int
  selectors::Vector{ColumnSelector}
  processes::Vector{GeoStatsProcess}
  rng::R
  kwargs::K
end

Simulate(domain::Domain, nreals::Int, selectors, processes, rng, kwargs) =
  Simulate(domain, nreals, collect(ColumnSelector, selectors), collect(GeoStatsProcess, processes), rng, kwargs)

Simulate(geoms::AbstractVector{<:Geometry}, nreals::Int, selectors, processes, rng, kwargs) =
  Simulate(GeometrySet(geoms), nreals, selectors, processes, rng, kwargs)

Simulate(domain, nreals::Int, pairs::Pair{<:Any,<:GeoStatsProcess}...; rng=Random.default_rng(), kwargs...) =
  Simulate(domain, nreals, selector.(first.(pairs)), last.(pairs), rng, values(kwargs))

Simulate(domain, pairs::Pair{<:Any,<:GeoStatsProcess}...; rng=Random.default_rng(), kwargs...) =
  Simulate(domain, 1, selector.(first.(pairs)), last.(pairs), rng, values(kwargs))

isrevertible(::Type{<:Simulate}) = false

function apply(transform::Simulate, geotable::AbstractGeoTable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  (; domain, nreals, selectors, processes, rng, kwargs) = transform
  ensembles = map(selectors, processes) do selector, process
    svars = selector(vars)
    data = geotable[:, svars]
    svars => rand(rng, process, domain, data, nreals; kwargs...)
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
