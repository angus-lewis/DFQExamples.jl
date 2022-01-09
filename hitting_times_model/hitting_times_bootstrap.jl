include((@__DIR__)*"/../preamble.jl") 

include("absorbing_model/model_def.jl")
include("default_params.jl")

for ic_string in ["exp","point_mass"]
    model_string = "qbdrap4"

    order = 1
    approx_data = CSV.read("hitting_times_model/hitting_times/data/"*
        ic_string*"/order_"*string(order)*"_model_"*model_string*".csv",DataFrame)

    sim_data = CSV.read("hitting_times_model/hitting_times/data/"*ic_string*"/sims.csv",DataFrame)

    function ht_cdf(sims)
        n_sims=length(sims.t)
        function F(x,i)
            i_idx = sims.φ.==i
            Fxi = sum(sims.t[i_idx].<=x)/n_sims
            return Fxi
        end
        return F
    end

    rng = StableRNGs.StableRNG(23122021)

    pth = mkpath("hitting_times_model/hitting_times/data/"*ic_string*"/bootstrap")
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_cdf_fun = ht_cdf(sim_data[idx,:])
        boot_cdf_evaluated = DataFrame(t=approx_data.t, 
                                        phase_1=boot_cdf_fun.(approx_data.t,1),
                                        phase_2=boot_cdf_fun.(approx_data.t,2))
        CSV.write(pth*"/bootsample_"*string(n)*".csv",boot_cdf_evaluated)
        print(string(n)*"... ")
    end
end


