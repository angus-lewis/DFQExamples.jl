include((@__DIR__)*"/../preamble.jl") 

include("absorbing_model/model_def.jl")
include("default_params.jl")

point_mass_initial_condition(n,rng) = (X = fill(eps(),n), φ = fill(1,n))
function exp_initial_condition(n,rng) 
    i = rand(rng,1:2,n) 
    X = Vector{Float64}(undef,n)
    for k in 1:n
        x = Inf
        while x>=1.0
            x = -log(rand(rng))
        end
        X[k]=x
    end
    return (X=X,φ=i)
end
rng = StableRNGs.StableRNG(16092021)

_first_exit = first_exit_x(0.5*eps(),1.0)
function first_exit_or_t10(model,SFM,SFM0)
    t1 = _first_exit(model,SFM,SFM0)
    i2 = SFM.t>10.0
    return (Ind=(i2||t1.Ind), SFM=t1.SFM)
end

sims = simulate(model,first_exit_or_t10,exp_initial_condition(n_sims,rng),rng)
pwd()
CSV.write("hitting_times_model/hitting_times/data/exp/sims.csv",
    DataFrame(t=sims.t,φ=sims.φ,X=sims.X,n=sims.n)    
)