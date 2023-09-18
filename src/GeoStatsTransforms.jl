# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsTransforms

using Meshes
using GeoTables
using GeoStatsModels

using Tables
using TableDistances
using TableTransforms
using ScientificTypes
using Combinatorics
using Distances
using Clustering
using ArnoldiMethod
using CategoricalArrays
using SparseArrays
using LinearAlgebra
using Statistics

using GeoStatsModels: GeoStatsModel, fit, predict, predictprob
import TableTransforms: ColSpec, Col, AllSpec, NoneSpec
import TableTransforms: colspec, choose
import TableTransforms: divide, attach
import TableTransforms: applymeta, revertmeta
import TableTransforms: apply, revert, reapply
import TableTransforms: isrevertible

include("utils.jl")

include("traits.jl")
include("feature.jl")
include("geometric.jl")

include("interpolate.jl")
include("uniquecoords.jl")
include("clustering.jl")
include("rasterize.jl")
include("potrace.jl")
include("detrend.jl")

export
  # transforms
  Interpolate,
  UniqueCoords,
  Rasterize,
  Potrace,
  Detrend,
  SLIC,
  GHC,
  GSC

end
