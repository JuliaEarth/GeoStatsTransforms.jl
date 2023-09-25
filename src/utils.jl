# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    searcher_ui(domain, maxneighbors, distance, neighborhood)

Return the appropriate search method over the `domain` based on
end-user inputs such as `maxneighbors`, `distance` and `neighborhood`.
"""
function searcher_ui(domain, maxneighbors, distance, neighborhood)
  if isnothing(neighborhood)
    # nearest neighbor search with a metric
    KNearestSearch(domain, maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(domain, maxneighbors, neighborhood)
  end
end

#-------------
# AGGREGATION
#-------------

defaultagg(x) = defaultagg(nonmissingtype(elscitype(x)), nonmissingtype(eltype(x)))
defaultagg(::Type{<:Continuous}, ::Type) = _mean
defaultagg(::Type, ::Type{<:AbstractQuantity}) = _mean
defaultagg(::Type, ::Type) = _first

function _mean(x)
  vs = skipmissing(x)
  isempty(vs) ? missing : mean(vs)
end

function _first(x)
  vs = skipmissing(x)
  isempty(vs) ? missing : first(vs)
end

#-------
# UNITS
#-------

function uadjust(geotable::AbstractGeoTable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  pairs = (var => uadjust(Tables.getcolumn(cols, var)) for var in vars)
  newtab = (; pairs...) |> Tables.materializer(tab)
  georef(newtab, dom)
end

uadjust(x) = uadjust(elunit(x), x)
uadjust(::Units, x) = x
function uadjust(u::AffineUnits, x)
  a = absoluteunit(u)
  [ismissing(v) ? missing : uconvert(a, v) for v in x]
end

elunit(x) = typeunit(nonmissingtype(eltype(x)))

typeunit(::Type) = NoUnits
typeunit(::Type{Q}) where {Q<:AbstractQuantity} = unit(Q)
