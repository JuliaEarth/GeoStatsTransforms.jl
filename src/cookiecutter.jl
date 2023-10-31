# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct CookieCutter{D<:Domain,M,C,R<:AbstractRNG,K} <: TableTransform
  domain::D
  nreals::Int
  master::M
  childrens::Vector{C}
  rng::R
  kwargs::K
end

function CookieCutter(
  domain::Domain,
  nreals::Int,
  master::Pair{C,<:GeoStatsProcess},
  childrens::Pair{C}...;
  rng=Random.default_rng(),
  kwargs...
) where {C<:Column}
  mpair = selector(first(master)) => last(master)
  cpairs = [selector(first(p)) => _childmap(last(p)) for p in childrens]
  CookieCutter(domain, nreals, mpair, cpairs, rng, values(kwargs))
end

function CookieCutter(
  domain::Domain,
  master::Pair{C,<:GeoStatsProcess},
  childrens::Pair{C}...;
  rng=Random.default_rng(),
  kwargs...
) where {C<:Column}
  mpair = selector(first(master)) => last(master)
  cpairs = [selector(first(p)) => _childmap(last(p)) for p in childrens]
  CookieCutter(domain, 1, mpair, cpairs, rng, values(kwargs))
end

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

  (; domain, nreals, master, childrens, rng, kwargs) = transform
  pad = ndigits(nreals)

  mselector, mprocess = master
  mvar = selectsingle(mselector, vars)
  mdata = geotable[:, mvar]
  mensemble = rand(rng, mprocess, domain, mdata, nreals; kwargs...)
  msim = mapreduce(hcat, 1:nreals) do r
    gtb = mensemble[r]
    gtb |> Rename(mvar => "$(mvar)_$(string(r; pad))")
  end

  csim = mapreduce(hcat, childrens) do (cselector, cmap)
    cvar = selectsingle(cselector, vars)
    cdata = geotable[:, cvar]
    prep = map(cmap) do (val, process)
      ensemble = rand(rng, process, domain, cdata, nreals; kwargs...)
      val => ensemble
    end
    censemble = Dict(prep)
    cpairs = map(1:nreals) do r
      mcolumn = msim[:, r]
      name = Symbol("$(cvar)_$(string(r; pad))")
      column = map(enumerate(mcolumn)) do (i, v)
        if haskey(censemble, v)
          ensemble = censemble[v]
          gtb = ensemble[r]
          gtb[:, cvar][i]
        else
          missing
        end
      end
      name => column
    end
    georef((; cpairs...), domain)
  end

  newgeotable = hcat(msim, csim)
  newgeotable, nothing
end
