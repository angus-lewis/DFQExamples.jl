import Pkg
println(pwd())
Pkg.activate(".")
if !occursin("DFQExamples.jl", Base.active_project())
    # Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
    # Pkg.rm("DiscretisedFluidQueues")
    display(!occursin("DFQExamples.jl", Base.active_project()))
    Pkg.develop(url=pwd()*"/../DiscretisedFluidQueues.jl")
    Pkg.instantiate()
end
using DiscretisedFluidQueues, Plots, DataFrames, CSV, LinearAlgebra, JLD2
import Distributions 

plotlyjs()