# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsTransforms

using CoordRefSystems
using Meshes
using GeoTables
using GeoStatsFunctions
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
using Random

using OhMyThreads: tmap
using TiledIteration: TileIterator
using ImageFiltering: KernelFactors, imgradients
using ColumnSelectors: ColumnSelector, SingleColumnSelector
using ColumnSelectors: Column, AllSelector, NoneSelector
using ColumnSelectors: selector, selectsingle
using GeoStatsModels: GeoStatsModel, fitpredict

import TableTransforms: apply, revert, reapply
import TableTransforms: isrevertible

include("utils.jl")

include("interpolate.jl")
include("interpneighbors.jl")
include("droplocallowhigh.jl")
include("uniquecoords.jl")
include("aggregate.jl")
include("transfer.jl")
include("upscale.jl")
include("downscale.jl")
include("gradient.jl")
include("rasterize.jl")
include("potrace.jl")
include("detrend.jl")
include("quenching.jl")
include("maxposterior.jl")
include("modefilter.jl")
include("clustering.jl")

include("precompile.jl")

export
  # transforms
  Interpolate,
  InterpolateNeighbors,
  DropLocalLowHigh,
  DropLocalLow,
  DropLocalHigh,
  UniqueCoords,
  Aggregate,
  Transfer,
  Upscale,
  Downscale,
  Gradient,
  Rasterize,
  Potrace,
  Detrend,
  Quenching,
  MaxPosterior,
  ModeFilter,
  SLIC,
  GHC,
  GSC

end
