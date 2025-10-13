# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Gradient(var; kern=KernelFactors.sobel)

Compute the gradient of given variable `var` over a grid.
Optionally, specify the `kern`el function from the package
[ImageFiltering.jl](https://github.com/JuliaImages/ImageFiltering.jl).

## Examples

```julia
Gradient(1) # gradient of first variable
Gradient("altitude") # slope of terrain
```
"""
struct Gradient{S<:ColumnSelector,K} <: TableTransform
  selector::S
  kern::K
end

Gradient(col::Column; kern=KernelFactors.sobel) = Gradient(selector(col), kern)

function apply(transform::Gradient, geotable::AbstractGeoTable)
  tab = values(geotable)
  dom = domain(geotable)

  # make sure parent domain is a regular grid
  dom isa RegularGrid || throw(ArgumentError("gradient only defined over regular grids"))

  # select target variable
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)
  svars = transform.selector(vars)

  # make sure only one variable is selected
  length(svars) == 1 || throw(ArgumentError("more than one variable selected"))

  # compute image gradient
  svar = first(svars)
  kern = transform.kern
  vals = Tables.getcolumn(cols, svar)
  grad = imgradients(reshape(vals, size(dom)), kern)

  # normalize by grid spacing
  ngrad = vec.(grad ./ spacing(dom))

  # define derivative names
  cnames = CoordRefSystems.names(crs(dom))
  newvar(i) = Symbol(svar, :_, cnames[i])
  newvars = ntuple(newvar, length(cnames))

  # construct attribute table
  newtable = (; zip(newvars, ngrad)...)

  # georeference over same domain
  newgeotable = georef(newtable, dom)

  newgeotable, nothing
end
