import Pkg
println(pwd())
Pkg.activate("../.")
# Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
# Pkg.rm("DiscretisedFluidQueues")
Pkg.develop(url=pwd()*"/../../DiscretisedFluidQueues.jl")
Pkg.instantiate(verbose=true)
using DiscretisedFluidQueues, Plots
# plotly()
gr()