using GeoStatsTransforms
using Meshes
using GeoTables
using Tables
using GeoStatsModels
using Variography
using TableTransforms
using ScientificTypes
using CategoricalArrays
using Statistics
using CoDa
using Test, Random
using FileIO: load

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  "feature.jl",
  "geometric.jl",
  "interpolate.jl",
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
