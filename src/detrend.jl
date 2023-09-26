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
  table = values(geotable)
  cols = Tables.columns(table)
  names = Tables.columnnames(cols)
  snames = transform.selector(names)

  tdata = trend(geotable, snames; degree=transform.degree)
  ttable = values(tdata)
  tcols = Tables.columns(ttable)

  ncols = map(names) do n
    x = Tables.getcolumn(cols, n)
    if n ∈ snames
      μ = Tables.getcolumn(tcols, n)
      x .- μ
    else
      x
    end
  end

  𝒯 = (; zip(names, ncols)...)
  newtable = 𝒯 |> Tables.materializer(table)

  newgeotable = georef(newtable, domain(geotable))

  newgeotable, (snames, tcols)
end

function revert(::Detrend, newgeotable, cache)
  newtable = values(newgeotable)
  cols = Tables.columns(newtable)
  names = Tables.schema(newtable).names

  snames, tcols = cache

  ncols = map(names) do n
    x = Tables.getcolumn(cols, n)
    if n ∈ snames
      μ = Tables.getcolumn(tcols, n)
      x .+ μ
    else
      x
    end
  end

  𝒯 = (; zip(names, ncols)...)
  table = 𝒯 |> Tables.materializer(newtable)

  georef(table, domain(newgeotable))
end

"""
    polymat(xs, d)

Return the matrix of monomials for the iterator `xs`, i.e.
for each item `x = (x₁, x₂,…, xₙ)` in `xs`, evaluate
the monomial terms of the expansion `(x₁ + x₂ + ⋯ + xₙ)ᵈ`
for a given degree `d`.

The resulting matrix has a number of rows that is equal
to the number of items in the iterator `xs`. The number
of columns is a function of the degree. For `d=0`, a
single column of ones is returned that corresponds to
the constant term `x₁⁰⋅x₂⁰⋅⋯⋅xₙ⁰` for all items in `xs`.
"""
function polymat(xs, d)
  x = first(xs)
  n = length(x)
  es = Iterators.flatten(multiexponents(n, d) for d in 0:d)
  ps = [[prod(x .^ e) for x in xs] for e in es]
  reduce(hcat, ps)
end

"""
    trend(data, vars; degree=1)

Return the deterministic spatial trend for the variables `vars`
in the spatial `data`. Approximate the trend with a polynomial
of given `degree`.

## References

* Menafoglio, A., Secchi, P. 2013. [A Universal Kriging predictor
  for spatially dependent functional data of a Hilbert Space]
  (https://doi.org/10.1214/13-EJS843)
"""
function trend(data, vars::AbstractVector{Symbol}; degree=1)
  𝒯 = values(data)
  𝒟 = domain(data)

  # retrieve columns
  cols = Tables.columns(𝒯)

  # build polynomial drift terms
  coords(𝒟, i) = coordinates(centroid(𝒟, i))
  xs = (coords(𝒟, i) for i in 1:nelements(𝒟))
  F = polymat(xs, degree)

  # eqs 25 and 26 in Menafoglio, A., Secchi, P. 2013.
  ms = map(vars) do var
    z = Tables.getcolumn(cols, var)
    a = (F'F \ F') * z
    F * a
  end

  ctor = Tables.materializer(𝒯)
  means = ctor((; zip(vars, ms)...))

  georef(means, 𝒟)
end

trend(data, var::Symbol; kwargs...) = trend(data, [var]; kwargs...)

function trend(data; kwargs...)
  𝒯 = values(data)
  s = Tables.schema(𝒯)
  vars = collect(s.names)
  trend(data, vars; kwargs...)
end
