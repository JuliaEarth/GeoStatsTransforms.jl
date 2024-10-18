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

using OhMyThreads: tmap
using Unitful: AffineQuantity
using TiledIteration: TileIterator
using ColumnSelectors: ColumnSelector, SingleColumnSelector
using ColumnSelectors: Column, AllSelector, NoneSelector
using ColumnSelectors: selector, selectsingle
using GeoStatsModels: GeoStatsModel, fitpredict
using GeoStatsProcesses: GeoStatsProcess

import TableTransforms: apply, revert, reapply
import TableTransforms: isrevertible

include("utils.jl")

include("interpolate.jl")
include("interpneighbors.jl")
include("interpmissing.jl")
include("interpnan.jl")
include("simulate.jl")
include("cookiecutter.jl")
include("uniquecoords.jl")
include("aggregate.jl")
include("transfer.jl")
include("upscale.jl")
include("downscale.jl")
include("clustering.jl")
include("rasterize.jl")
include("potrace.jl")
include("detrend.jl")

include("precompile.jl")

export
  # transforms
  InterpolateNeighbors,
  Interpolate,
  InterpolateMissing,
  InterpolateNaN,
  Simulate,
  CookieCutter,
  UniqueCoords,
  Aggregate,
  Transfer,
  Upscale,
  Downscale,
  Rasterize,
  Potrace,
  Detrend,
  SLIC,
  GHC,
  GSC

end
