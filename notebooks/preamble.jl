import Pkg
Pkg.activate("..")
Pkg.instantiate()
Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
using DiscretisedFluidQueues, Plots
pyplot()