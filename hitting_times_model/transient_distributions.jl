# ## Run this from the hitting_times_model directory
# using Distributed
# if nprocs()==1
#     addprocs(4)
# end
# @everywhere begin 
#     include("everywhere_transient_distributions.jl")
# end

# # run the jobs 
# @sync begin 
#     t1 = @spawnat :any run_transient("reflecting_model","point_mass",point_mass_initial_condition,d0_map_point_mass,StableRNGs.StableRNG(16092021))

#     t2 = @spawnat :any run_transient("reflecting_model","exp",exp_initial_condition,d0_map_exp,StableRNGs.StableRNG(16092021))

#     t3 = @spawnat :any run_transient("absorbing_model","point_mass",point_mass_initial_condition,d0_map_point_mass,StableRNGs.StableRNG(16092021))

#     t4 = @spawnat :any run_transient("absorbing_model","exp",exp_initial_condition,d0_map_exp,StableRNGs.StableRNG(16092021))
# end 


## comment out the above, and uncomment the below for a non-paralell test 
include("everywhere_transient_distributions.jl")
run_transient("absorbing_model","point_mass",exp_initial_condition,d0_map_exp,StableRNGs.StableRNG(16092021))