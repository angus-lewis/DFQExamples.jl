# include("../preamble.jl")
# include((@__DIR__)*"/../../preamble.jl")

T = [-1.1 1.1; 1.0 -1.0]
c = [1.0; -1.0]
b_absorbing = 10.0 

T_absorbing = [T [0.0 0.0; 0.0 0.0]; [0.01 0.0 -0.01 0.0; 0.0 0.01 0.0 -0.01]]
c_absorbing = [c; 0.0; 0.0]

P_upr_absorbing = [0.0 0.0 0.0 1.0]
P_lwr_absorbing = [0.0 0.0 1.0 0.0] # both absorbing boundaries

model = BoundedFluidQueue(T_absorbing,c_absorbing,P_lwr_absorbing,P_upr_absorbing,b_absorbing)


