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

_addunit(x::Number, u) = x * u
_addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))

_adjustunits(geotable::AbstractGeoTable) = georef(_adjustunits(values(geotable)), domain(geotable))

function _adjustunits(table)
  cols = Tables.columns(table)
  vars = Tables.columnnames(cols)

  pairs = (var => _absunit(Tables.getcolumn(cols, var)) for var in vars)
  (; pairs...) |> Tables.materializer(table)
end

_absunit(x) = _absunit(nonmissingtype(eltype(x)), x)
_absunit(::Type, x) = x
function _absunit(::Type{Q}, x) where {Q<:AffineQuantity}
  u = absoluteunit(unit(Q))
  map(v -> uconvert(u, v), x)
end
