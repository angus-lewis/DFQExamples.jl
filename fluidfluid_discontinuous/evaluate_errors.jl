include("default_params.jl")
include("model_def.jl")

ks_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_cell_probs_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
ks_error_row = zeros(4)
l1_error_row = zeros(4)
l1_cell_probs_row = zeros(4)

sim_cdf = CSV.read(
    (@__DIR__)*"/data/sims_cdf_evaluated.csv",
    DataFrame,
)

for o in orders
    cols = [1;3;4]
    for (c,m) in enumerate(string.(approx_types[[1;3;4]]).*string.([1;3;4]))
        approx_data = CSV.read(
            (@__DIR__)*"/data/first_return_cdf_approximations/order_"*string(o)*"_model_"*m*".csv",
            DataFrame,
        )
        ks_error = max(
            maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
            maximum(abs.(approx_data.phase_4-sim_cdf.phase_4)),
        )
        ks_error_row[cols[c]] = ks_error
        l1_cdf = sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.x))/2 + 
            sum(abs.(approx_data.phase_4[1:end-1]+approx_data.phase_4[2:end]-
                sim_cdf.phase_4[1:end-1]-sim_cdf.phase_4[2:end]).*diff(approx_data.x))/2
        l1_error_row[cols[c]] = l1_cdf
        cell_edges_idx = findall(mod.(approx_data.x,0.4).==0.0)
        sim_cell_probs = DataFrame(
            x=[approx_data.x[1]; approx_data.x[cell_edges_idx[2:end]] - approx_data.x[cell_edges_idx[1:end-1]]],
            phase_2=[sim_cdf.phase_2[1]; sim_cdf.phase_2[cell_edges_idx[2:end]] - sim_cdf.phase_2[cell_edges_idx[1:end-1]]],
            phase_4=[sim_cdf.phase_4[1]; sim_cdf.phase_4[cell_edges_idx[2:end]] - sim_cdf.phase_4[cell_edges_idx[1:end-1]]],
        )
        approx_cell_probs = DataFrame(
            x=[approx_data.x[1]; approx_data.x[cell_edges_idx[2:end]] - approx_data.x[cell_edges_idx[1:end-1]]],
            phase_2=[approx_data.phase_2[1]; approx_data.phase_2[cell_edges_idx[2:end]] - approx_data.phase_2[cell_edges_idx[1:end-1]]],
            phase_4=[approx_data.phase_4[1]; approx_data.phase_4[cell_edges_idx[2:end]] - approx_data.phase_4[cell_edges_idx[1:end-1]]],
        )
        l1_cell_probs_row[cols[c]] = sum(abs.(approx_cell_probs.phase_2-sim_cell_probs.phase_2)) + 
            sum(abs.(approx_cell_probs.phase_4-sim_cell_probs.phase_4))
    end
    push!(ks_errors,ks_error_row)
    push!(l1_errors,l1_error_row)
    push!(l1_cell_probs_errors,l1_cell_probs_row)
end

CSV.write((@__DIR__)*"/data/errors/ks_error.csv",ks_errors)
CSV.write((@__DIR__)*"/data/errors/l1_error.csv",l1_errors)
CSV.write((@__DIR__)*"/data/errors/l1_cell_probs_error.csv",l1_cell_probs_errors)

ks_errors_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_errors_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_cell_probs_errors_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
ks_error_row_lwr = zeros(4)
l1_error_row_lwr = zeros(4)
l1_cell_probs_error_row_lwr = zeros(4)

ks_errors_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_errors_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_cell_probs_errors_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
ks_error_row_upr = zeros(4)
l1_error_row_upr = zeros(4)
l1_cell_probs_error_row_upr = zeros(4)

ks_error_boot = zeros(n_boot)
l1_error_boot = zeros(n_boot)
l1_cell_probs_error_boot = zeros(n_boot)
for o in orders
    cols = [1;3;4]
    for (c,m) in enumerate(string.(approx_types[[1;3;4]]).*string.([1;3;4]))
        approx_data = CSV.read(
            (@__DIR__)*"/data/first_return_cdf_approximations/order_"*string(o)*"_model_"*m*".csv",
            DataFrame,
        )
        cell_edges_idx = findall(mod.(approx_data.x,0.4).==0.0)
        approx_cell_probs = DataFrame(
                x=[approx_data.x[1]; approx_data.x[cell_edges_idx[2:end]] - approx_data.x[cell_edges_idx[1:end-1]]],
                phase_2=[approx_data.phase_2[1]; approx_data.phase_2[cell_edges_idx[2:end]] - approx_data.phase_2[cell_edges_idx[1:end-1]]],
                phase_4=[approx_data.phase_4[1]; approx_data.phase_4[cell_edges_idx[2:end]] - approx_data.phase_4[cell_edges_idx[1:end-1]]],
            )
        for n in 1:n_boot
            sim_cdf = CSV.read(
                (@__DIR__)*"/data/bootstrap_first_return/sample_"*string(n)*".csv",
                DataFrame,
            )
            ks_error = max(
                maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
                maximum(abs.(approx_data.phase_4-sim_cdf.phase_4)),
            )
            ks_error_boot[n] = ks_error
            l1_cdf = sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                    sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.x))/2 + 
                sum(abs.(approx_data.phase_4[1:end-1]+approx_data.phase_4[2:end]-
                    sim_cdf.phase_4[1:end-1]-sim_cdf.phase_4[2:end]).*diff(approx_data.x))/2
            l1_error_boot[n] = l1_cdf
            sim_cell_probs = DataFrame(
                x=[approx_data.x[1]; approx_data.x[cell_edges_idx[2:end]] - approx_data.x[cell_edges_idx[1:end-1]]],
                phase_2=[sim_cdf.phase_2[1]; sim_cdf.phase_2[cell_edges_idx[2:end]] - sim_cdf.phase_2[cell_edges_idx[1:end-1]]],
                phase_4=[sim_cdf.phase_4[1]; sim_cdf.phase_4[cell_edges_idx[2:end]] - sim_cdf.phase_4[cell_edges_idx[1:end-1]]],
            )
            l1_cell_probs_error_boot[n] = sum(abs.(approx_cell_probs.phase_2-sim_cell_probs.phase_2)) + 
                sum(abs.(approx_cell_probs.phase_4-sim_cell_probs.phase_4))
        end
        ks_error_row_lwr[cols[c]] = quantile(ks_error_boot,qtiles[1])
        l1_error_row_lwr[cols[c]] = quantile(l1_error_boot,qtiles[1])
        l1_cell_probs_error_row_lwr[cols[c]] = quantile(l1_cell_probs_error_boot,qtiles[1])
        ks_error_row_upr[cols[c]] = quantile(ks_error_boot,qtiles[2])
        l1_error_row_upr[cols[c]] = quantile(l1_error_boot,qtiles[2])
        l1_cell_probs_error_row_upr[cols[c]] = quantile(l1_cell_probs_error_boot,qtiles[2])
    end
    push!(ks_errors_lwr,ks_error_row_lwr)
    push!(l1_errors_lwr,l1_error_row_lwr)
    push!(l1_cell_probs_errors_lwr,l1_cell_probs_error_row_lwr)
    push!(ks_errors_upr,ks_error_row_upr)
    push!(l1_errors_upr,l1_error_row_upr)
    push!(l1_cell_probs_errors_upr,l1_cell_probs_error_row_upr)
end

CSV.write((@__DIR__)*"/data/errors/ks_error_lwr.csv",ks_errors_lwr)
CSV.write((@__DIR__)*"/data/errors/l1_error_lwr.csv",l1_errors_lwr)
CSV.write((@__DIR__)*"/data/errors/l1_cell_probs_error_lwr.csv",l1_cell_probs_errors_lwr)
CSV.write((@__DIR__)*"/data/errors/ks_error_upr.csv",ks_errors_upr)
CSV.write((@__DIR__)*"/data/errors/l1_error_upr.csv",l1_errors_upr)
CSV.write((@__DIR__)*"/data/errors/l1_cell_probs_error_upr.csv",l1_cell_probs_errors_upr)

