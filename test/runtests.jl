using GeoStatsTransforms
using Meshes
using Tables
using Unitful
using GeoTables
using GeoStatsFunctions
using GeoStatsModels
using GeoStatsProcesses
using GeoStatsImages
using TableTransforms
using CategoricalArrays
using Statistics
using Test, StableRNGs
using FileIO: load
import DataScienceTraits as DST

import ImageQuilting

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  "interpolate.jl",
  "interpneighbors.jl",
  "interpmissing.jl",
  "interpnan.jl",
  "simulate.jl",
  "cookiecutter.jl",
  "uniquecoords.jl",
  "aggregate.jl",
  "transfer.jl",
  "upscale.jl",
  "downscale.jl",
  "clustering.jl",
  "rasterize.jl",
  "potrace.jl",
  "detrend.jl"
]

@testset "GeoStatsTransforms.jl" begin
  for testfile in testfiles
    println("Testing $testfile...")
    include(testfile)
  end
end
