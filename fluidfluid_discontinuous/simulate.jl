include("default_params.jl")
include("model_def.jl")

seed = 10012022
rng = StableRNGs.StableRNG(seed)
# reset!(rng) = Random.seed!(rng,seed)

# n_sims=493079
sfm_sims, y_sims = simulate(
    ffq_sim, # model to sim
    first_exit_y(0.0,ffq_sim.dq.model.b+1.0), # stopping time (first return to 0.0)
    (φ = fill(4,n_sims), X = fill(2.0,n_sims), Y = zeros(n_sims)), # initial condition
    rng,
)

sim_data = DataFrame((sfm_sims...,Y=y_sims))

CSV.write((@__DIR__)*"/data/sims.csv",sim_data)
