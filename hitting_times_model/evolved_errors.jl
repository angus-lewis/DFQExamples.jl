include((@__DIR__)*"/../preamble.jl") 

include("reflecting_model/model_def.jl")
include("default_params.jl")

pth = mkpath("hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/errors")

ic_string = "point_mass"

let 
    idx = [1
    1001
    2001
    3001
    4001
    5001
    6001
    7001
    8001
    9001
   10001]
    ks_error = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    l1_error = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    cell_probs_error = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    ks_error_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    l1_error_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    cell_probs_error_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    ks_error_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    l1_error_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
    cell_probs_error_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])

    sim_cdf = CSV.read("hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/"*
                    "sims_cdf_evaluated.csv",DataFrame)
    for order in orders
        ks_error_row = zeros(4)
        l1_error_row = zeros(4)
        cell_probs_error_row = zeros(4)
        ks_error_row_lwr = zeros(4)
        l1_error_row_lwr = zeros(4)
        cell_probs_error_row_lwr = zeros(4)
        ks_error_row_upr = zeros(4)
        l1_error_row_upr = zeros(4)
        cell_probs_error_row_upr = zeros(4)
        for (c,model_string) in enumerate(string.([approx_types...]).*string.(1:4))
            approx_data = CSV.read("hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/approx_cdf/order_"*string(order)*"model_"*model_string*".csv",DataFrame)

            ks_error_row[c] = max(
                maximum(abs.(approx_data.phase_1-sim_cdf.phase_1)),
                maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
            )
            l1_error_row[c] = sum(abs.(approx_data.phase_1[1:end-1]+approx_data.phase_1[2:end]-
                    sim_cdf.phase_1[1:end-1]-sim_cdf.phase_1[2:end]).*diff(approx_data.x))/2 + 
                sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                    sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.x))/2
            
            cell_probs_error_row[c] = sum(abs.(diff(approx_data.phase_1[idx])-diff(sim_cdf.phase_1[idx]))) +
                sum(abs.(diff(approx_data.phase_2[idx])-diff(sim_cdf.phase_2[idx])))

            ks_error_samples = zeros(n_boot)
            l1_error_samples = zeros(n_boot)
            cell_probs_error_samples = zeros(n_boot)
            for n in 1:n_boot
                cdf_boot = CSV.read("hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/bootstrap/sample_"*string(n)*".csv",DataFrame)
                
                ks_error_samples[n] = max(
                    maximum(abs.(approx_data.phase_1-cdf_boot.phase_1)),
                    maximum(abs.(approx_data.phase_2-cdf_boot.phase_2)),
                )
                l1_error_samples[n] = sum(abs.(approx_data.phase_1[1:end-1]+approx_data.phase_1[2:end]-
                        cdf_boot.phase_1[1:end-1]-cdf_boot.phase_1[2:end]).*diff(approx_data.x))/2 + 
                    sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                        cdf_boot.phase_2[1:end-1]-cdf_boot.phase_2[2:end]).*diff(approx_data.x))/2
                
                cell_probs_error_samples[n] = sum(abs.(diff(approx_data.phase_1[idx])-diff(cdf_boot.phase_1[idx]))) +
                    sum(abs.(diff(approx_data.phase_2[idx])-diff(cdf_boot.phase_2[idx])))
            end
            ks_q5 = quantile(ks_error_samples,qtiles[1])
            ks_q95 = quantile(ks_error_samples,qtiles[2])
            l1_q5 = quantile(l1_error_samples,qtiles[1])
            l1_q95 = quantile(l1_error_samples,qtiles[2])
            cp_q5 = quantile(cell_probs_error_samples,qtiles[1])
            cp_q95 = quantile(cell_probs_error_samples,qtiles[2])
            ks_error_row_lwr[c] = ks_q5
            ks_error_row_upr[c] = ks_q95
            l1_error_row_lwr[c] = l1_q5
            l1_error_row_upr[c] = l1_q95
            cell_probs_error_row_lwr[c] = cp_q5
            cell_probs_error_row_upr[c] = cp_q95
        end
        push!(ks_error_lwr,ks_error_row_lwr)
        push!(ks_error_upr,ks_error_row_upr)
        push!(l1_error_lwr,l1_error_row_lwr)
        push!(l1_error_upr,l1_error_row_upr)
        push!(cell_probs_error_lwr,cell_probs_error_row_lwr)
        push!(cell_probs_error_upr,cell_probs_error_row_upr)
        push!(ks_error,ks_error_row)
        push!(l1_error,l1_error_row)
        push!(cell_probs_error,cell_probs_error_row)
    end
    CSV.write(pth*"/ks_lwr_"*ic_string*".csv",ks_error_lwr)
    CSV.write(pth*"/ks_upr_"*ic_string*".csv",ks_error_upr)
    CSV.write(pth*"/l1_lwr_"*ic_string*".csv",l1_error_lwr)
    CSV.write(pth*"/l1_upr_"*ic_string*".csv",l1_error_upr)
    CSV.write(pth*"/cell_probs_lwr_"*ic_string*".csv",cell_probs_error_lwr)
    CSV.write(pth*"/cell_probs_upr_"*ic_string*".csv",cell_probs_error_upr)
    CSV.write(pth*"/ks_"*ic_string*".csv",ks_error_lwr)
    CSV.write(pth*"/l1_"*ic_string*".csv",l1_error_lwr)
    CSV.write(pth*"/cell_probs_"*ic_string*".csv",cell_probs_error)
end
