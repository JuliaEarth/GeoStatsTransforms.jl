using GeoStatsTransforms
using Meshes
using Tables
using Unitful
using GeoTables
using Variography
using GeoStatsModels
using GeoStatsProcesses
using TableTransforms
using CategoricalArrays
using Statistics
using Test, Random
using FileIO: load
import DataScienceTraits as DST

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  "feature.jl",
  "geometric.jl",
  "interpneighbors.jl",
  "interpolate.jl",
  "simulate.jl",
  "uniquecoords.jl",
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
