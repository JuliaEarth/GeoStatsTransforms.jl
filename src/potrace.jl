# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Potrace(mask; [ϵ])
    Potrace(mask, var₁ => agg₁, ..., varₙ => aggₙ; [ϵ])

Trace polygons on 2D image data with Selinger's Potrace algorithm.

The categories stored in column `mask` are converted into binary
masks, which are then traced into multi-polygons. When provided,
the option `ϵ` is forwarded to Selinger's simplification algorithm.

Duplicates of a variable `varᵢ` are aggregated with
aggregation function `aggᵢ`. If an aggregation function 
is not defined for variable `varᵢ`, the default aggregation 
function will be used. Default aggregation function is `mean` for
continuous variables and `first` otherwise.

## Examples

```julia
Potrace(:mask, ϵ=0.1)
Potrace(1, 1 => last, 2 => maximum)
Potrace(:mask, :a => first, :b => minimum)
Potrace("mask", "a" => last, "b" => maximum)
```

## References

- Selinger, P. 2003. [Potrace: A polygon-based tracing algorithm]
  (https://potrace.sourceforge.net/potrace.pdf)
"""
struct Potrace{M<:SingleColumnSelector,S<:ColumnSelector,T} <: TableTransform
  mask::M
  selector::S
  aggfuns::Vector{Function}
  ϵ::T
end

Potrace(mask::Column; ϵ=nothing) = Potrace(selector(mask), NoneSelector(), Function[], ϵ)
Potrace(mask::Column, pairs::Pair{C,<:Function}...; ϵ=nothing) where {C<:Column} =
  Potrace(selector(mask), selector(first.(pairs)), collect(Function, last.(pairs)), ϵ)

isrevertible(::Type{<:Potrace}) = true

function apply(transform::Potrace, geotable::AbstractGeoTable)
  gtb = geotable |> AbsoluteUnits()
  dom = domain(gtb)
  tab = values(gtb)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  # sanity check
  if !(dom isa Grid)
    throw(ArgumentError("potrace only defined for grid data"))
  end

  # simplification threshold
  ϵ = transform.ϵ

  # select column name
  sname = selectsingle(transform.mask, vars)

  # aggregation functions
  svars = transform.selector(vars)
  agg = Dict(zip(svars, transform.aggfuns))
  for var in vars
    if !haskey(agg, var)
      v = Tables.getcolumn(cols, var)
      agg[var] = _defaultagg(v)
    end
  end

  # convert column to image
  col = Tables.getcolumn(cols, sname)
  img = reshape(col, size(dom))

  # all possible colors
  colors = unique(img)

  # aggregate variables within each color
  preproc = map(colors) do color
    mask = isequal.(img, color)
    inds = findall(vec(mask))
    feat = Any[sname => color]
    for var in setdiff(vars, [sname])
      v = Tables.getcolumn(cols, var)
      newv = agg[var](v[inds])
      push!(feat, var => newv)
    end
    (; feat...), mask
  end

  # split preprocessing results
  feats = first.(preproc)
  masks = last.(preproc)

  # collect vertices and topology
  verts = vertices(dom)
  topo = topology(dom)

  # map pixels to vertices
  ∂ = Boundary{2,0}(topo)

  # map direction to first vertex of edge of
  # interest (i.e. edge touched by direction)
  d = Dict(:→ => 1, :↑ => 2, :← => 3, :↓ => 4)

  # map (→, i) representation to ring of points
  ring(itr) = Ring([verts[∂(i)[d[→]]] for (→, i) in itr[begin:(end - 1)]])

  # trace multi-polygons on each mask
  multis = map(masks) do mask
    rings = trace(mask)
    polys = map(rings) do (outer, inners)
      ochain = ring(outer)
      ichains = [ring(inner) for inner in inners]
      PolyArea([ochain, ichains...])
    end
    Multi(polys)
  end

  # simplify multi-polygons if necessary
  elems = isnothing(ϵ) ? multis : [simplify(multi, SelingerSimplification(ϵ)) for multi in multis]

  # georeference new features on new geometries
  newtab = feats |> Tables.materializer(tab)
  newdom = elems |> GeometrySet
  newgeotable = georef(newtab, newdom)

  newgeotable, dom
end

revert(::Potrace, newgeotable::AbstractGeoTable, cache) = newgeotable |> Rasterize(cache)

# trace polygonal geometries on mask
function trace(mask)
  # pad mask with inactive pixels
  M = falses(size(mask) .+ 2)
  M[(begin + 1):(end - 1), (begin + 1):(end - 1)] .= mask

  # trace paths on padded mask
  paths = tracerecursion!(M)

  # convert paths into rings
  rings = paths2rings(paths)

  # unpad and linearize indices
  linear = LinearIndices(mask)
  fun(■) = linear[■ - CartesianIndex(1, 1)]
  map(rings) do (outer, inners)
    o = [(→, fun(■)) for (□, →, ■) in outer]
    is = [[(→, fun(■)) for (□, →, ■) in inner] for inner in inners]
    o, is
  end
end

function tracerecursion!(M)
  paths = []
  while any(M)
    # trace outer path
    outer = tracepath(M)

    # invert pixels inside path
    O = copy(M)
    insideout!(M, outer)
    I = @. M & !O
    @. M = M & !I

    if any(I)
      # perform recursion
      inners = tracerecursion!(I)
      push!(paths, (outer, inners))
    else
      # single outer path
      push!(paths, (outer, []))
    end
  end

  paths
end

# trace the top-left polygon on the mask
function tracepath(M)
  # find top-left corner (□ → ■ link)
  i, j = 1, findfirst(==(1), M[1, :])
  while isnothing(j) && i < size(M, 1)
    i += 1
    j = findfirst(==(1), M[i, :])
  end

  # there must be at least one active pixel
  @assert !isnothing(j) "invalid input mask"

  # define □ → ■ link
  □ = CartesianIndex(i, j - 1)
  ■ = CartesianIndex(i, j)

  # step direction along the path
  step(□, ■) = CartesianIndex(■[2] - □[2], □[1] - ■[1])

  # direction after a given turn
  left = Dict(:→ => :↑, :↑ => :←, :← => :↓, :↓ => :→)
  right = Dict(:→ => :↓, :↓ => :←, :← => :↑, :↑ => :→)

  # find the next edge along the path
  function move((□, →, ■))
    □ₛ = □ + step(□, ■)
    ■ₛ = ■ + step(□, ■)

    # 4 possible configurations
    if M[□ₛ] == 1 && M[■ₛ] == 1
      □, right[→], □ₛ # make a right turn
    elseif M[□ₛ] == 0 && M[■ₛ] == 1
      □ₛ, →, ■ₛ # continue straight
    elseif M[□ₛ] == 0 && M[■ₛ] == 0
      ■ₛ, left[→], ■ # make a left turn
    else # cross pattern
      ■ₛ, left[→], ■ # left turn policy
    end
  end

  # build a closed path
  start = (□, :→, ■)
  next = move(start)
  path = [start, next]
  while next ≠ start
    next = move(next)
    push!(path, next)
  end

  path
end

# invert the the mask inside the path
function insideout!(M, path)
  □s, ⬕s = first.(path), last.(path)
  frontier = collect(zip(□s, ⬕s))
  visited = falses(size(M))
  visited[□s] .= true
  while !isempty(frontier)
    □, ⬕ = pop!(frontier)

    if !visited[⬕]
      # flip color
      M[⬕] = 1 - M[⬕]
      visited[⬕] = true

      # update frontier
      δ = ⬕ - □
      ⬕₁ = ⬕ + δ
      ⬕₂ = ⬕ + CartesianIndex(δ[2], -δ[1])
      ⬕₃ = ⬕ + CartesianIndex(-δ[2], δ[1])
      for ⬕ₛ in [⬕₁, ⬕₂, ⬕₃]
        if !visited[⬕ₛ]
          push!(frontier, (⬕, ⬕ₛ))
        end
      end
    end
  end

  M
end

# convert forest of paths to rings
paths2rings(paths) = mapreduce(treebfs, vcat, paths)

# breadth-first-search on tree of paths
function treebfs(root)
  record(node) = (first(node), first.(last(node)))

  visited = []
  frontier = [root]
  while !isempty(frontier)
    node = popfirst!(frontier)

    seen = false
    for vnode in visited
      if first(node) ∈ last(vnode)
        seen = true
        break
      end
    end

    if !seen
      push!(visited, record(node))
    end

    for child in last(node)
      push!(frontier, child)
    end
  end

  visited
end
