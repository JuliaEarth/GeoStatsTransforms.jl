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

  # parent domain and indices
  grid = parent(dom)
  inds = parentindices(dom)

  # make sure parent domain is a regular grid
  grid isa RegularGrid || throw(ArgumentError("gradient only defined over views of regular grids"))

  # retrieve grid spacing
  spac = spacing(grid)

  # select target variable
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)
  svars = transform.selector(vars)

  # make sure only one variable is selected
  length(svars) == 1 || throw(ArgumentError("more than one variable selected for gradient calculation"))

  # reshape column into image format
  svar = first(svars)
  vals = Tables.getcolumn(cols, svar)
  imag = zeros(eltype(vals), size(grid))
  imag[inds] .= vals

  # compute image gradient
  kern = transform.kern
  grad = imgradients(imag, kern)

  # normalize by grid spacing
  ngrad = map(zip(grad, spac)) do (∂, h)
    ∂[inds] / h
  end

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
