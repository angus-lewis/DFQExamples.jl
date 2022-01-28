using CSV
using DataFrames
using Plots 
using JLD2
using DiscretisedFluidQueues
import Distributions
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/error_metrics.jl")

@halfwidth_plot_defaults()

plotlyjs()

ks_data = CSV.read((@__DIR__)*"/../errors/ks_point_mass.csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../errors/ks_lwr_point_mass.csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../errors/ks_upr_point_mass.csv",DataFrame)
@show lwr = [log10(ks_data[i,j]) .- log10(ks_data_lwr[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
@show upr = [log10(ks_data_upr[i,j]) .- log10(ks_data[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
ks_data = log10.(ks_data)
p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:false,:false,:false,:false]
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
colours = [1;3;2;4]
for (c,col) in enumerate(names(ks_data))
    plot!(1:2:21,ks_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="KS error")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/ks_error_formatted.pdf")

@halfwidth_plot_defaults()

L1_cell_probs_data = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_point_mass.csv",DataFrame))
L1_cell_probs_data_lwr = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_lwr_point_mass.csv",DataFrame))
L1_cell_probs_data_upr = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_upr_point_mass.csv",DataFrame))
@show lwr = [(L1_cell_probs_data[i,j]) .- (L1_cell_probs_data_lwr[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
@show upr = [(L1_cell_probs_data_upr[i,j]) .- (L1_cell_probs_data[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]


p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:false,:false,:false,:false]
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
colours = [1;3;2;4]
for (c,col) in enumerate(names(L1_cell_probs_data))
    plot!(1:2:21,L1_cell_probs_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="Cell-wise error")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/L1_cell_probs_error_formatted.pdf")



l1_cdf_data = log10.(CSV.read((@__DIR__)*"/../errors/l1_point_mass.csv",DataFrame))
l1_cdf_data_lwr = log10.(CSV.read((@__DIR__)*"/../errors/l1_lwr_point_mass.csv",DataFrame))
l1_cdf_data_upr = log10.(CSV.read((@__DIR__)*"/../errors/l1_upr_point_mass.csv",DataFrame))

@show lwr = [(l1_cdf_data[i,j]) .- (l1_cdf_data_lwr[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
@show upr = [(l1_cdf_data_upr[i,j]) .- (l1_cdf_data[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]

@halfwidth_plot_defaults()
p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
for (c,col) in enumerate(names(l1_cdf_data))
    plot!(1:2:21,l1_cdf_data[:,col],ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="L¹ error - CDF")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")

