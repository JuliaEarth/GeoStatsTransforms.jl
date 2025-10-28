# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DropLocalLowHigh(radius; low=0.25, high=0.75)

Drop rows of geotable when values are outside percentile
interval `[low, high]` within a given `radius`.

## Examples

```julia
# drop rows outside interquartile range using a 10m radius
DropLocalLowHigh(10u"m")

# drop rows with extremely large values using 5m radius
DropLocalLowHigh(5u"m", low=0.0, high=0.98)
```

See also [`DropLocalMinima`](@ref) and [`DropLocalMaxima`](@ref).
"""
struct DropLocalLowHigh{ℒ<:Len,T} <: TableTransform
  radius::ℒ
  low::T
  high::T
end

DropLocalLowHigh(radius; low=0.25, high=0.75) = DropLocalLowHigh(aslen(radius), low, high)

function apply(transform::DropLocalLowHigh, geotable::AbstractGeoTable)
  # domain and continuous variables
  dom = domain(geotable)
  tab = values(geotable) |> Only(Continuous)

  # columns and variable names
  cols = Tables.columns(tab)
  vars = Tables.columnnames(tab)

  # transform parameters
  radius = transform.radius
  low = transform.low
  high = transform.high

  # define search method
  ball = MetricBall(radius)
  searcher = BallSearch(dom, ball)

  # find rows to drop
  drop = Int[]
  for i in 1:nelements(dom)
    inds = search(centroid(dom, i), searcher)
    for var in vars
      x = Tables.getcolumn(cols, var)
      l, h = quantile(x[inds], (low, high))
      if x[i] < l || x[i] > h
        push!(drop, i)
        break
      end
    end
  end

  # find rows to keep
  keep = setdiff(1:nelements(dom), drop)

  # return subgeotable
  geotable[keep, :], nothing
end

"""
    DropLocalMinima(radius; low=0.25)

Equivalent to `DropLocalLowHigh(radius; low=low, high=1.0)`.

Please check the documentation of [`DropLocalLowHigh`](@ref) for more details.
"""
DropLocalMinima(radius; low=0.25) = DropLocalLowHigh(radius; low=low, high=1.0)

"""
    DropLocalMaxima(radius; high=0.75)

Equivalent to `DropLocalLowHigh(radius; low=0.0, high=high)`.

Please check the documentation of [`DropLocalLowHigh`](@ref) for more details.
"""
DropLocalMaxima(radius; high=0.75) = DropLocalLowHigh(radius; low=0.0, high=high)
