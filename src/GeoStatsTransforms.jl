# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsTransforms

using Meshes
using GeoTables
using GeoStatsModels
using GeoStatsProcesses

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
using Random

using Unitful: AffineQuantity
using GeoStatsModels: GeoStatsModel, fitpredict
using ColumnSelectors: ColumnSelector, SingleColumnSelector
using ColumnSelectors: Column, AllSelector, NoneSelector
using ColumnSelectors: selector, selectsingle
using GeoStatsProcesses: GeoStatsProcess
using DataScienceTraits: Continuous

import TableTransforms: apply, revert, reapply
import TableTransforms: isrevertible

include("utils.jl")

include("interpneighbors.jl")
include("interpolate.jl")
include("simulate.jl")
include("cookiecutter.jl")
include("uniquecoords.jl")
include("aggregate.jl")
include("transfer.jl")
include("upscale.jl")
include("clustering.jl")
include("rasterize.jl")
include("potrace.jl")
include("detrend.jl")

export
  # transforms
  InterpolateNeighbors,
  Interpolate,
  Simulate,
  CookieCutter,
  UniqueCoords,
  Aggregate,
  Transfer,
  Upscale,
  Rasterize,
  Potrace,
  Detrend,
  SLIC,
  GHC,
  GSC

end
