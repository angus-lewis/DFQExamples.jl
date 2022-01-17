include("default_params.jl")
include("model_def.jl")

sim_data = CSV.read((@__DIR__)*"/data/sims.csv",DataFrame)
x_vals=range(0,1.6;length=n_err_evals)

function cdf(x::AbstractVector,data::DataFrame) 
    sort!(data,:X)
    
    x2 = data.X[data.φ.==2]
    x4 = data.X[data.φ.==4]

    c2 = 1
    c4 = 1

    F2 = zeros(length(x))
    F4 = zeros(length(x))

    for (c,xi) in enumerate(x)
        tmp2 = findfirst(y->y>xi,x2[c2:end])
        if !(tmp2===nothing)
            c2 += tmp2-1
        else 
            c2 = length(x2)+1
        end
        F2[c] = c2-1
        
        tmp4 = findfirst(y->y>xi,x4[c4:end])
        if !(tmp4===nothing)
            c4 += tmp4-1
        else 
            c4 = length(x4)+1
        end
        F4[c] = c4-1
    end
    return F2./length(data.X), F4./length(data.X)
end

F2,F4 = cdf(x_vals,sim_data)

CSV.write((@__DIR__)*"/data/sims_cdf_evaluated.csv",DataFrame(phase_2=F2,phase_4=F4))

function bootstrap(sim_data,x_vals,n_boot)
    rng = StableRNGs.StableRNG(1412)
    n_sims = length(sim_data.X)
    for n in 1:n_boot
        ind = rand(rng,1:n_sims,n_sims)
        boot_data = sim_data[ind,:]
        F2,F4 = cdf(x_vals,boot_data)
        CSV.write(
            (@__DIR__)*"/data/bootstrap_first_return/sample_"*string(n)*".csv",
            DataFrame(phase_2=F2,phase_4=F4),
        )
    end
    return nothing 
end

bootstrap(sim_data,x_vals,n_boot)