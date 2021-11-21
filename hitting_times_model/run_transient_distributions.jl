## Run this from the hitting_times_model directory
using Distributed
# if nprocs()==1
addprocs(4)
# end
@everywhere begin 
    include("transient_distributions.jl")
end

# run the jobs 
@sync begin 
    t1 = @spawnat :any begin
        include("reflecting_model/model_def.jl")
        include("default_params.jl")
        run_transient(
            orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
            model,"reflecting_model",point_mass_initial_condition,"point_mass",d0_map_point_mass,StableRNGs.StableRNG(16092021)
        )
    end

    t2 = @spawnat :any begin
        include("reflecting_model/model_def.jl")
        include("default_params.jl")
        run_transient(
            orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
            model,"reflecting_model",exp_initial_condition,"exp",d0_map_exp,StableRNGs.StableRNG(16092021)
        )
    end

    t3 = @spawnat :any begin 
        include("absorbing_model/model_def.jl")
        include("default_params.jl")
        run_transient(
            orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
            model,"absorbing_model",point_mass_initial_condition,"point_mass",d0_map_point_mass,StableRNGs.StableRNG(16092021)
        )
    end 

    t4 = @spawnat :any begin 
        include("absorbing_model/model_def.jl")
        include("default_params.jl")
        run_transient(
            orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
            model,"absorbing_model",exp_initial_condition,"exp",d0_map_exp,StableRNGs.StableRNG(16092021)
        )
    end
end 


## comment out the above, and uncomment the below for a non-paralell test 
# include("transient_distributions.jl")
# include("absorbing_model/model_def.jl")
# include("default_params.jl")
# run_transient(
#     orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
#     model,"absorbing_model",exp_initial_condition,"exp",d0_map_exp,StableRNGs.StableRNG(16092021)
# )