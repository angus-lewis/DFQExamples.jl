import Pkg
println(@__DIR__)
if !occursin("DFQExamples.jl", Base.active_project())
    # Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
    # Pkg.rm("DiscretisedFluidQueues")
    Pkg.activate(@__DIR__)
    display(!occursin("DFQExamples.jl", Base.active_project()))
    Pkg.develop(url=(@__DIR__)*"/../DiscretisedFluidQueues.jl")
    Pkg.instantiate()
end
using DiscretisedFluidQueues, Plots, DataFrames, CSV, LinearAlgebra, JLD2
import Distributions 
include((@__DIR__)*"/SFM_operators.jl")

plotlyjs()