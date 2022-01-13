include("default_params.jl")
include("model_def.jl")

ks_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
l1_errors = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
ks_error_row = zeros(4)
l1_error_row = zeros(4)

sim_cdf = CSV.read(
    (@__DIR__)*"/data/sims_cdf_evaluated.csv",
    DataFrame,
)

for o in orders 
    for (c,m) in enumerate(string.(approx_types).*string.(1:4))
        approx_data = CSV.read(
            (@__DIR__)*"/data/first_return_cdf_approximations/order_"*string(o)*"_model_"*m*".csv",
            DataFrame,
        )
        ks_error = max(
            maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
            maximum(abs.(approx_data.phase_4-sim_cdf.phase_4)),
        )
        ks_error_row[c] = ks_error
        l1_cdf = sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
                sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.x))/2 + 
            sum(abs.(approx_data.phase_4[1:end-1]+approx_data.phase_4[2:end]-
                sim_cdf.phase_4[1:end-1]-sim_cdf.phase_4[2:end]).*diff(approx_data.x))/2
        l1_error_row[c] = l1_cdf
    end
    push!(ks_errors,ks_error_row)
    push!(l1_errors,l1_error_row)
end

CSV.write((@__DIR__)*"/data/errors/ks_error.csv",ks_errors)
CSV.write((@__DIR__)*"/data/errors/l1_error.csv",l1_errors)

# ks_errors_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
# l1_errors_lwr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
# ks_error_row_lwr = zeros(4)
# l1_error_row_lwr = zeros(4)

# ks_errors_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
# l1_errors_upr = DataFrame(DG=[],DG_limit=[],Order_1=[],QBDRAP=[])
# ks_error_row_upr = zeros(4)
# l1_error_row_upr = zeros(4)

# ks_error_boot = zeros(n_boot)
# l1_error_boot = zeros(n_boot)
# for o in orders 
#     for (c,m) in enumerate(string.(approx_types).*string.(1:4))
#         approx_data = CSV.read(
#             (@__DIR__)*"/data/first_return_cdf_approximations/order_"*string(o)*"_model_"*m*".csv",
#             DataFrame,
#         )
#         for n in 1:n_boot
#             sim_cdf = CSV.read(
#                 (@__DIR__)*"/data/bootstrap_first_return/sample_"*string(n)*".csv",
#                 DataFrame,
#             )
#             ks_error = max(
#                 maximum(abs.(approx_data.phase_2-sim_cdf.phase_2)),
#                 maximum(abs.(approx_data.phase_4-sim_cdf.phase_4)),
#             )
#             ks_error_boot[n] = ks_error
#             l1_cdf = sum(abs.(approx_data.phase_2[1:end-1]+approx_data.phase_2[2:end]-
#                     sim_cdf.phase_2[1:end-1]-sim_cdf.phase_2[2:end]).*diff(approx_data.x))/2 + 
#                 sum(abs.(approx_data.phase_4[1:end-1]+approx_data.phase_4[2:end]-
#                     sim_cdf.phase_4[1:end-1]-sim_cdf.phase_4[2:end]).*diff(approx_data.x))/2
#             l1_error_boot[n] = l1_cdf
#         end
#         ks_error_row_lwr[c] = quantile(ks_error_boot,qtiles[1])
#         l1_error_row_lwr[c] = quantile(l1_error_boot,qtiles[1])
#         ks_error_row_upr[c] = quantile(ks_error_boot,qtiles[2])
#         l1_error_row_upr[c] = quantile(l1_error_boot,qtiles[2])
#     end
#     push!(ks_errors_lwr,ks_error_row_lwr)
#     push!(l1_errors_lwr,l1_error_row_lwr)
#     push!(ks_errors_upr,ks_error_row_upr)
#     push!(l1_errors_upr,l1_error_row_upr)
# end

# CSV.write((@__DIR__)*"/data/errors/ks_error_lwr.csv",ks_errors_lwr)
# CSV.write((@__DIR__)*"/data/errors/l1_error_lwr.csv",l1_errors_lwr)
# CSV.write((@__DIR__)*"/data/errors/ks_error_upr.csv",ks_errors_upr)
# CSV.write((@__DIR__)*"/data/errors/l1_error_upr.csv",l1_errors_upr)