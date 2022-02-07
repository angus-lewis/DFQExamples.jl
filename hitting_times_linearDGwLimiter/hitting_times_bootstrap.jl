include((@__DIR__)*"/../preamble.jl") 

include("absorbing_model/model_def.jl")
include("default_params.jl")

for ic_string in ["point_mass"]
    model_string = "qbdrap4"

    order = 1
    approx_data = CSV.read("hitting_times_model/hitting_times/data/"*
        ic_string*"/order_"*string(order)*"_model_"*model_string*".csv",DataFrame)
    x_vals = approx_data.t

    sim_data = CSV.read("hitting_times_model/hitting_times/data/"*ic_string*"/sims.csv",DataFrame)

    function cdf(x::AbstractVector,data::DataFrame) 
        sort!(data,:t)
        
        x1 = data.t[data.φ.==1]
        x2 = data.t[data.φ.==2]
    
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
        return F1./length(data.t), F2./length(data.t)
    end
    
    F1,F2 = cdf(x_vals,sim_data)
    
    CSV.write((@__DIR__)*"/hitting_times/data/"*ic_string*"/sims_cdf_evaluated.csv",DataFrame(t=x_vals,phase_1=F1,phase_2=F2))
    
    function bootstrap(sim_data,x_vals,n_boot)
        rng = StableRNGs.StableRNG(1412)
        n_sims = length(sim_data.t)
        for n in 1:n_boot
            @show n
            ind = rand(rng,1:n_sims,n_sims)
            boot_data = sim_data[ind,:]
            F1,F2 = cdf(x_vals,boot_data)
            CSV.write(
                (@__DIR__)*"/hitting_times/data/"*ic_string*"/bootstrap/bootsample_"*string(n)*".csv",
                DataFrame(t=x_vals,phase_1=F1,phase_2=F2),
            )
        end
        return nothing 
    end
    
    bootstrap(sim_data,x_vals,n_boot)

    # rng = StableRNGs.StableRNG(23122021)

    # pth = mkpath("hitting_times_model/hitting_times/data/"*ic_string*"/bootstrap")
    # for n in 1:n_boot
    #     idx = rand(rng,1:n_sims,n_sims)
    #     boot_cdf_fun = ht_cdf(sim_data[idx,:])
    #     boot_cdf_evaluated = DataFrame(t=approx_data.t, 
    #                                     phase_1=boot_cdf_fun.(approx_data.t,1),
    #                                     phase_2=boot_cdf_fun.(approx_data.t,2))
    #     CSV.write(pth*"/bootsample_"*string(n)*".csv",boot_cdf_evaluated)
    #     print(string(n)*"... ")
    # end
end


