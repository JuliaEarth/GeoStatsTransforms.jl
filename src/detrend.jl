# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Detrend(col₁, col₂, ..., colₙ; degree=1)
    Detrend([col₁, col₂, ..., colₙ]; degree=1)
    Detrend((col₁, col₂, ..., colₙ); degree=1)
    
The transform that detrends columns `col₁`, `col₂`, ..., `colₙ`
with a polynomial of given `degree`.

    Detrend(regex; degree=1)

Detrends the columns that match with `regex`.

# Examples

```julia
Detrend(1, 3, 5)
Detrend([:a, :c, :e])
Detrend(("a", "c", "e"))
Detrend(r"[ace]", degree=2)
Detrend(:)
```

## References

* Menafoglio, A., Secchi, P. 2013. [A Universal Kriging predictor
  for spatially dependent functional data of a Hilbert Space]
  (https://doi.org/10.1214/13-EJS843)
"""
struct Detrend{S<:ColumnSelector} <: TableTransform
  selector::S
  degree::Int
end

Detrend(; degree=1) = Detrend(AllSelector(), degree)
Detrend(cols; degree=1) = Detrend(selector(cols), degree)
Detrend(cols::C...; degree=1) where {C<:Column} = Detrend(selector(cols), degree)

isrevertible(::Type{<:Detrend}) = true

function apply(transform::Detrend, geotable)
  dom = domain(geotable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  names = Tables.columnnames(cols)
  snames = transform.selector(names)

  gview = geotable |> Select(snames)
  model = Polynomial(transform.degree)
  fitted = GeoStatsModels.fit(model, gview)

  ncols = map(names) do name
    z = Tables.getcolumn(cols, name)
    ẑ(i) = GeoStatsModels.predict(fitted, name, centroid(dom, i))
    if name ∈ snames
      @inbounds [z[i] - ẑ(i) for i in 1:nelements(dom)]
    else
      z
    end
  end

  𝒯 = (; zip(names, ncols)...)
  newtab = 𝒯 |> Tables.materializer(tab)

  newgeotable = georef(newtab, dom)

  newgeotable, (snames, fitted)
end

function revert(::Detrend, newgeotable, cache)
  newdom = domain(newgeotable)
  newtab = values(newgeotable)
  cols = Tables.columns(newtab)
  names = Tables.columnnames(cols)

  snames, fitted = cache

  ocols = map(names) do name
    z = Tables.getcolumn(cols, name)
    ẑ(i) = GeoStatsModels.predict(fitted, name, centroid(newdom, i))
    if name ∈ snames
      @inbounds [z[i] + ẑ(i) for i in 1:nelements(newdom)]
    else
      z
    end
  end

  𝒯 = (; zip(names, ocols)...)
  table = 𝒯 |> Tables.materializer(newtab)

  georef(table, newdom)
end
