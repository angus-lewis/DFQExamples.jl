include("default_params.jl")
include("model_def.jl")

seed = 10012022
rng = StableRNGs.StableRNG(seed)
# reset!(rng) = Random.seed!(rng,seed)

n_sims=20_000
sfm_sims, y_sims = simulate(
    ffq_sim, # model to sim
    first_exit_y(0.0,ffq_sim.dq.model.b+1.0), # stopping time (first return to 0.0)
    (φ = fill(3,n_sims), X = fill(5.0,n_sims), Y = zeros(n_sims)), # initial condition
    rng,
)

sim_data = DataFrame((sfm_sims...,Y=y_sims))

# CSV.write((@__DIR__)*"/data/sims.csv",sim_data)

# f(x) = 0.01.+x^2
# FluidFluidQueues.brute_force_zero(f,-1,1)

# x = range(-1,1;step=1e-6)
# fx = abs.(f.(x))
