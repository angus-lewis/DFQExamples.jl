include((@__DIR__)*"/../preamble.jl") 

include("reflecting_model/model_def.jl")
include("default_params.jl")

pth = mkpath((@__DIR__)*"/hitting_times/data/errors")

for ic_string in (ic,)
    ks_error_lwr = DataFrame(DG_limit=[])
    l1_error_lwr = DataFrame(DG_limit=[])
    ks_error_upr = DataFrame(DG_limit=[])
    l1_error_upr = DataFrame(DG_limit=[])
    for order in orders
        ks_error_row_lwr = zeros(1)
        l1_error_row_lwr = zeros(1)
        ks_error_row_upr = zeros(1)
        l1_error_row_upr = zeros(1)
        # for (c,model_string) in [:dg;]
        c=1
            ks_error_samples = zeros(n_boot)
            l1_error_samples = zeros(n_boot)
            for n in 1:n_boot
                cdf_boot = CSV.read("hitting_times_model/hitting_times/data/"*
                    ic_string*"/bootstrap/bootsample_"*string(n)*".csv",DataFrame)

                approx_data = CSV.read((@__DIR__)*"/hitting_times/data/"*
                    ic_string*"/order_"*string(order)*"_model_dg1.csv",DataFrame)
                
                ks_error_samples[n] = max(
                    maximum(abs.(approx_data.phase_1-cdf_boot.phase_1)),
                    maximum(abs.(approx_data.phase_2-cdf_boot.phase_2)),
                )
                l1_error_samples[n] = sum(abs.(approx_data.phase_1[1:end-1]+approx_data.phase_1[2:end]-
                        cdf_boot.phase_1[1:end-1]-cdf_boot.phase_1[2:end]).*diff(approx_data.t))/2 + 
                    sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                        cdf_boot.phase_2[1:end-1]-cdf_boot.phase_2[2:end]).*diff(approx_data.t))/2
            end
            ks_q5 = quantile(ks_error_samples,qtiles[1])
            ks_q95 = quantile(ks_error_samples,qtiles[2])
            l1_q5 = quantile(l1_error_samples,qtiles[1])
            l1_q95 = quantile(l1_error_samples,qtiles[2])
            ks_error_row_lwr[c] = ks_q5
            ks_error_row_upr[c] = ks_q95
            l1_error_row_lwr[c] = l1_q5
            l1_error_row_upr[c] = l1_q95
        # end
        push!(ks_error_lwr,ks_error_row_lwr)
        push!(ks_error_upr,ks_error_row_upr)
        push!(l1_error_lwr,l1_error_row_lwr)
        push!(l1_error_upr,l1_error_row_upr)
    end
    CSV.write(pth*"/ks_lwr_"*ic_string*".csv",ks_error_lwr)
    CSV.write(pth*"/ks_upr_"*ic_string*".csv",ks_error_upr)
    CSV.write(pth*"/l1_lwr_"*ic_string*".csv",l1_error_lwr)
    CSV.write(pth*"/l1_upr_"*ic_string*".csv",l1_error_upr)
end
