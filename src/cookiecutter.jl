# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CookieCutter(domain, master => process, var₁ => map₁, ..., varₙ => mapₙ; [parameters])
    CookieCutter(domain, nreals, master => process, var₁ => map₁, ..., varₙ => mapₙ; [parameters])

TODO

## Parameters

* `rng`      - Random number generator (default to `Random.default_rng()`)
* `pool`     - Pool of worker processes (default to `[myid()]`)
* `threads`  - Number of threads (default to `cpucores()`)
* `progress` - Show progress bar (default to `true`)
"""
struct CookieCutter{D<:Domain,M,C,R<:AbstractRNG,K} <: TableTransform
  domain::D
  nreals::Int
  master::M
  children::Vector{C}
  rng::R
  kwargs::K
end

function CookieCutter(
  domain::Domain,
  nreals::Int,
  master::Pair{C,<:GeoStatsProcess},
  children::Pair{C}...;
  rng=Random.default_rng(),
  kwargs...
) where {C<:Column}
  if isempty(children)
    throw(ArgumentError("cannot create CookieCutter transform without children"))
  end
  mpair = selector(first(master)) => last(master)
  cpairs = [selector(first(p)) => _childmap(last(p)) for p in children]
  CookieCutter(domain, nreals, mpair, cpairs, rng, values(kwargs))
end

CookieCutter(domain::Domain, master::Pair{C,<:GeoStatsProcess}, children::Pair{C}...; kwargs...) where {C<:Column} =
  CookieCutter(domain, 1, master, children...; kwargs...)

function _childmap(itr)
  pairs = collect(itr)
  if !(eltype(pairs) <: Pair{<:Any,<:GeoStatsProcess})
    throw(ArgumentError("child map must be an iterable of pairs of the form: value => process"))
  end
  pairs
end

isrevertible(::Type{<:CookieCutter}) = false

function apply(transform::CookieCutter, geotable::AbstractGeoTable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  (; domain, nreals, master, children, rng, kwargs) = transform

  mselector, mprocess = master
  mvar = selectsingle(mselector, vars)
  msim = geotable |> Simulate(domain, nreals, mvar => mprocess; rng, kwargs...)

  csim = mapreduce(hcat, children) do (cselector, cpairs)
    cvar = selectsingle(cselector, vars)
    prep = map(cpairs) do (val, cprocess)
      sim = geotable |> Simulate(domain, nreals, cvar => cprocess; rng, kwargs...)
      val => sim
    end

    names = let
      tab = values(last(first(prep)))
      cols = Tables.columns(tab)
      Tables.columnnames(cols)
    end

    cmap = Dict(prep)
    columns = map(1:nreals) do r
      mcolumn = msim[:, r]
      map(enumerate(mcolumn)) do (i, v)
        if haskey(cmap, v)
          sim = cmap[v]
          sim[:, r][i]
        else
          missing
        end
      end
    end

    georef((; zip(names, columns)...), domain)
  end

  newgeotable = hcat(msim, csim)
  newgeotable, nothing
end
