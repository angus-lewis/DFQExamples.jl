include((@__DIR__)*"/../preamble.jl") 

include("absorbing_model/model_def.jl")
include("default_params.jl")

function ht_cdf(sims)
    n_sims=length(sims.t)
    function F(x,i)
        i_idx = sims.φ.==i
        Fxi = sum(sims.t[i_idx].<=x)/n_sims
        return Fxi
    end
    return F
end

tvec = CSV.read("hitting_times_model/hitting_times/data/point_mass/order_1_model_dg1.csv",DataFrame).t

pth = mkpath("hitting_times_model/hitting_times/data/errors")
for ic_string in ["point_mass","exp"]
    ks_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    l1_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    ks_error_row = zeros(4)
    l1_error_row = zeros(4)

    sim_data = CSV.read("hitting_times_model/hitting_times/data/"*ic_string*"/sims.csv",DataFrame)

    sim_cdf = DataFrame(
        t=tvec,
        phase_1=ht_cdf(sim_data).(tvec,1),
        phase_2=ht_cdf(sim_data).(tvec,2),
    )
    for order in 1:2:21
        for (c,model_string) in enumerate(["dg1","dg2","order13","qbdrap4"])
            approx_data = CSV.read("hitting_times_model/hitting_times/data/"*
                ic_string*"/order_"*string(order)*"_model_"*model_string*".csv",DataFrame)
            ks_error = max(
                maximum(abs.(approx_data.phase_1-sim_cdf.phase_1)),
                maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
            )
            l1_cdf = sum(abs.(approx_data.phase_1[1:end-1]+approx_data.phase_1[2:end]-
                    sim_cdf.phase_1[1:end-1]-sim_cdf.phase_1[2:end]).*diff(approx_data.t))/2 + 
                sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                    sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.t))/2
            
            ks_error_row[c] = ks_error
            l1_error_row[c] = l1_cdf
        end
        push!(ks_errors,ks_error_row)
        push!(l1_errors,l1_error_row)
    end
    CSV.write(pth*"/ks_"*ic_string*".csv",ks_errors)
    CSV.write(pth*"/l1_"*ic_string*".csv",l1_errors)
end