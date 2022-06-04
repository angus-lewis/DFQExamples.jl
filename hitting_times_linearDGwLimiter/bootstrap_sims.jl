include("default_params.jl")
include("reflecting_model/model_def.jl")

sim_data = CSV.read(
    (@__DIR__)*
    "/../hitting_times_model/reflecting_model/transient_distribution/"*
    ic*"/data/sim.csv",DataFrame)
x_vals=range(0,model.b;length=n_err_evals)

function cdf(x::AbstractVector,data::DataFrame) 
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

F1,F2 = cdf(x_vals,sim_data)

CSV.write((@__DIR__)*"/reflecting_model/data/"*ic*"/data/sims_cdf_evaluated.csv",DataFrame(phase_1=F1,phase_2=F2))

function bootstrap(sim_data,x_vals,n_boot)
    rng = StableRNGs.StableRNG(14129)
    n_sims = length(sim_data.X)
    for n in 1:n_boot
        ind = rand(rng,1:n_sims,n_sims)
        boot_data = sim_data[ind,:]
        F1,F2 = cdf(x_vals,boot_data)
        CSV.write(
            (@__DIR__)*"/reflecting_model/data/"*ic*"/data/bootstrap/sample_"*string(n)*".csv",
            DataFrame(phase_1=F1,phase_2=F2),
        )
    end
    return nothing 
end

bootstrap(sim_data,x_vals,n_boot)