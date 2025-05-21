using Graphs
using Random
using SimpleWeightedGraphs
using SparseArrays
using StatsBase

# missing other dependencies
# temporary
@info ("Temprary import structure... Build CommunityDetectionAlgorithms module")
include(joinpath("lib", "types.jl"))
include(joinpath("lib", "fixed_stack.jl"))
include(joinpath("lib", "lib.jl"))
include(joinpath("lib", "partitions.jl"))
include(joinpath("algorithms", "Leiden.jl"))

#=
module CommunityDetectionAlgorithms




end # module CommunityDetectionAlgorithms
=#