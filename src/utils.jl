# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#-------------
# AGGREGATION
#-------------

_defaultagg(x) = _defaultagg(elscitype(x))
_defaultagg(::Type) = _skipmissing(first)
_defaultagg(::Type{Continuous}) = _skipmissing(mean)

function _skipmissing(fun)
  x -> begin
    vs = skipmissing(x)
    isempty(vs) ? missing : fun(vs)
  end
end

#-------
# UNITS
#-------

const Len{T} = Quantity{T,u"𝐋"}

_addunit(x::Number, u) = x * u
_addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))

function _adjustunits(geotable::AbstractGeoTable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  pairs = (var => _absunit(Tables.getcolumn(cols, var)) for var in vars)
  newtab = (; pairs...) |> Tables.materializer(tab)

  georef(newtab, dom)
end

_absunit(x) = _absunit(nonmissingtype(eltype(x)), x)
_absunit(::Type, x) = x
function _absunit(::Type{Q}, x) where {Q<:AffineQuantity}
  u = absoluteunit(unit(Q))
  map(v -> uconvert(u, v), x)
end

#---------
# THREADS
#---------

function _tmap(f, itr)
  nchunks = cld(length(itr), Threads.nthreads())
  chunks = Iterators.partition(itr, nchunks)
  tasks = map(chunks) do chunk
    Threads.@spawn map(f, chunk)
  end
  mapreduce(fetch, vcat, tasks)
end
