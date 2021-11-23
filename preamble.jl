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

using CSV 
using DiscretisedFluidQueues 
using DataFrames 
using Distributed
import Distributions 
using Plots 
using LinearAlgebra 
using JLD2 
using JSON 
using Random 
using StableRNGs 
using Statistics


include((@__DIR__)*"/error_metrics.jl")
include((@__DIR__)*"/helper_functions.jl")

plotlyjs()