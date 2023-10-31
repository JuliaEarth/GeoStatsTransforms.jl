# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CookieCutter(domain, parent => process, var₁ => procmap₁, ..., varₙ => procmapₙ; [parameters])
    CookieCutter(domain, nreals, parent => process, var₁ => procmap₁, ..., varₙ => procmapₙ; [parameters])

Simulate `nreals` realizations of variable `parent` with geostatistical process `process`, 
and each child variable `varsᵢ` with process map `procmapᵢ`, over given `domain`.

The process map must be an iterable of pairs of the form: value => process.
Each process in the map is related to a value of the `parent` realization, 
therefore the values of the child variables will be chosen according
to the values of the corresponding `parent` realization.

## Parameters

* `rng`      - Random number generator (default to `Random.default_rng()`)
* `pool`     - Pool of worker processes (default to `[myid()]`)
* `threads`  - Number of threads (default to `cpucores()`)
* `progress` - Show progress bar (default to `true`)

# Examples

```julia
parent = QuiltingProcess(trainimg=trainimg, tilesize=(30, 30))
child0 = GaussianProcess(variogram=SphericalVariogram(range=20.0, sill=0.2))
child1 = GaussianProcess(variogram=SphericalVariogram(MetricBall((200.0, 20.0))))
transform = CookieCutter(domain, :parent_var => parent, :child_var => [0 => child0, 1 => child1])
```
"""
struct CookieCutter{D<:Domain,M,C,R<:AbstractRNG,K} <: TableTransform
  domain::D
  nreals::Int
  parent::M
  children::Vector{C}
  rng::R
  kwargs::K
end

function CookieCutter(
  domain::Domain,
  nreals::Int,
  parent::Pair{C,<:GeoStatsProcess},
  children::Pair{C}...;
  rng=Random.default_rng(),
  kwargs...
) where {C<:Column}
  if isempty(children)
    throw(ArgumentError("cannot create CookieCutter transform without children"))
  end
  ppair = selector(first(parent)) => last(parent)
  cpairs = [selector(first(p)) => _procmap(last(p)) for p in children]
  CookieCutter(domain, nreals, ppair, cpairs, rng, values(kwargs))
end

CookieCutter(domain::Domain, parent::Pair{C,<:GeoStatsProcess}, children::Pair{C}...; kwargs...) where {C<:Column} =
  CookieCutter(domain, 1, parent, children...; kwargs...)

function _procmap(itr)
  pairs = collect(itr)
  if !(eltype(pairs) <: Pair{<:Any,<:GeoStatsProcess})
    throw(ArgumentError("process map must be an iterable of pairs of the form: value => process"))
  end
  pairs
end

isrevertible(::Type{<:CookieCutter}) = false

function apply(transform::CookieCutter, geotable::AbstractGeoTable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  (; domain, nreals, parent, children, rng, kwargs) = transform

  pselector, pprocess = parent
  pvar = selectsingle(pselector, vars)
  psim = geotable |> Simulate(domain, nreals, pvar => pprocess; rng, kwargs...)

  csim = mapreduce(hcat, children) do (cselector, procmap)
    cvar = selectsingle(cselector, vars)
    prep = map(procmap) do (val, cprocess)
      sim = geotable |> Simulate(domain, nreals, cvar => cprocess; rng, kwargs...)
      val => sim
    end

    names = let
      tab = values(last(first(prep)))
      cols = Tables.columns(tab)
      Tables.columnnames(cols)
    end

    simmap = Dict(prep)
    columns = map(1:nreals) do r
      mcolumn = psim[:, r]
      map(enumerate(mcolumn)) do (i, v)
        if haskey(simmap, v)
          sim = simmap[v]
          sim[:, r][i]
        else
          missing
        end
      end
    end

    georef((; zip(names, columns)...), domain)
  end

  newgeotable = hcat(psim, csim)
  newgeotable, nothing
end
