import Pkg
println(pwd())
Pkg.activate("../.")
Pkg.instantiate(verbose=true)
# Pkg.add(url="https://github.com/angus-lewis/DiscretisedFluidQueues")
Pkg.develop(url=pwd()*"/../../DiscretisedFluidQueues")
using DiscretisedFluidQueues, Plots
# pyplot()