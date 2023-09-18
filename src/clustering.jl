# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClusteringTransform

A transform for geostatistical clustering.
"""
abstract type ClusteringTransform <: TableTransform end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clustering/slic.jl")
include("clustering/ghc.jl")
include("clustering/gsc.jl")
