defaultagg(x) = defaultagg(nonmissingtype(elscitype(x)))
defaultagg(::Type{<:Continuous}) = _mean
defaultagg(::Type) = _first

function _mean(x)
  vs = skipmissing(x)
  isempty(vs) ? missing : mean(vs)
end

function _first(x)
  vs = skipmissing(x)
  isempty(vs) ? missing : first(vs)
end

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
