import Pkg
println(@__DIR__)
if !occursin("DFQExamples.jl", Base.active_project())
    # Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
    # Pkg.rm("DiscretisedFluidQueues")
    Pkg.activate(@__DIR__)
    Pkg.develop(url=(@__DIR__)*"/../DiscretisedFluidQueues.jl")
    Pkg.develop(url=(@__DIR__)*"/../FluidQueues.jl")
    Pkg.instantiate()
    display(!occursin("DFQExamples.jl", Base.active_project()))
end
using DiscretisedFluidQueues 
using Plots 
using DataFrames 
using CSV 
using LinearAlgebra 
using JLD2 
using Random 
using StableRNGs 
using Statistics
using Distributed
import Distributions 

include((@__DIR__)*"/error_metrics.jl")
include((@__DIR__)*"/helper_functions.jl")

plotlyjs()