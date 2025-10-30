# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# fit tuple `dims` to a given length `D` by repeating the last dimension.
_fitdims(dims::Dims{N}, D) where {N} = ntuple(i -> i ≤ N ? dims[i] : last(dims), D)

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

aslen(x::Len) = x
aslen(x::Number) = x * u"m"
aslen(::Quantity) = throw(ArgumentError("invalid length unit"))

# ------
# STATS
# ------

function _mode(levs, vals)
  c = Dict(levs .=> 0)
  @inbounds for v in vals
    c[v] += 1
  end
  argmax(c)
end
