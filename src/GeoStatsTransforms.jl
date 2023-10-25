# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsTransforms

using Meshes
using GeoTables
using GeoStatsModels

using Tables
using Unitful
using TableDistances
using TableTransforms
using DataScienceTraits
using Combinatorics
using Distances
using Clustering
using ArnoldiMethod
using CategoricalArrays
using SparseArrays
using LinearAlgebra
using Statistics

using Unitful: AffineQuantity
using GeoStatsModels: GeoStatsModel, fitpredict
using ColumnSelectors: ColumnSelector, SingleColumnSelector
using ColumnSelectors: Column, AllSelector, NoneSelector
using ColumnSelectors: selector, selectsingle
using DataScienceTraits: Continuous

import TableTransforms: divide, attach
import TableTransforms: applymeta, revertmeta
import TableTransforms: apply, revert, reapply
import TableTransforms: isrevertible

include("utils.jl")

include("traits.jl")
include("feature.jl")
include("geometric.jl")

include("interpneighbors.jl")
include("interpolate.jl")
include("uniquecoords.jl")
include("clustering.jl")
include("rasterize.jl")
include("potrace.jl")
include("detrend.jl")

export
  # transforms
  InterpolateNeighbors,
  Interpolate,
  UniqueCoords,
  Rasterize,
  Potrace,
  Detrend,
  SLIC,
  GHC,
  GSC

end
