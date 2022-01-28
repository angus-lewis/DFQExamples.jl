include((@__DIR__)*"/../preamble.jl") 
include("reflecting_model/model_def.jl")

sims = CSV.read((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/data/sim.csv",DataFrame)

rng = StableRNGs.StableRNG(26092021)

sims = simulate(model,fixed_time(0.1),(X=sims.X,φ=sims.φ),rng)
pwd()
sim_data = DataFrame(t=sims.t,φ=sims.φ,X=sims.X,n=sims.n)  
CSV.write((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2/sims.csv",
    sim_data,    
)

x_vals=range(0,10;length=n_err_evals)

function fast_cdf(x::AbstractVector,data::DataFrame) 
    sort!(data,:X)
    
    x1 = data.X[data.φ.==1]
    x2 = data.X[data.φ.==2]

    c1 = 1
    c2 = 1

    F1 = zeros(length(x))
    F2 = zeros(length(x))

    for (c,xi) in enumerate(x)
        tmp1 = findfirst(y->y>xi,x1[c1:end])
        if !(tmp1===nothing)
            c1 += tmp1-1
        else 
            c1 = length(x1)+1
        end
        F1[c] = c1-1
        
        tmp2 = findfirst(y->y>xi,x2[c2:end])
        if !(tmp2===nothing)
            c2 += tmp2-1
        else 
            c2 = length(x2)+1
        end
        F2[c] = c2-1
    end
    return F1./length(data.X), F2./length(data.X)
end

F1, F2 = fast_cdf(x_vals,sim_data)

CSV.write((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2/sims_cdf_evaluated.csv",DataFrame(x=x_vals,phase_1=F1,phase_2=F2))

function bootstrap(sim_data,x_vals,n_boot)
    rng = StableRNGs.StableRNG(1412)
    n_sims = length(sim_data.X)
    for n in 1:n_boot
        @show n
        ind = rand(rng,1:n_sims,n_sims)
        boot_data = sim_data[ind,:]
        F1, F2 = fast_cdf(x_vals,boot_data)
        CSV.write(
            (@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2/bootstrap/sample_"*string(n)*".csv",
            DataFrame(x=x_vals,phase_1=F1,phase_2=F2),
        )
    end
    return nothing 
end

bootstrap(sim_data,x_vals,n_boot)