using CSV
using DataFrames
using Plots 
using JLD2
using DiscretisedFluidQueues
using GLM
import Distributions
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/error_metrics.jl")

@halfwidth_plot_defaults()

plotlyjs()
function p()
ks_data = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnov.csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnov_ci_lwr.csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnovci_upr.csv",DataFrame)
 lwr = [log10(ks_data[i,j]) .- log10(ks_data_lwr[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
 upr = [log10(ks_data_upr[i,j]) .- log10(ks_data[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
ks_data = log10.(ks_data)
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:false,:false,:false,:false]
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
colours = [1;3;2;4]
for (c,col) in enumerate(names(ks_data))
    plot!(log10.(1:2:21),ks_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error")
plot!(legend=:outerbottomright)
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/ks_error_formatted.pdf")
end 
p()

@halfwidth_plot_defaults()
function p()
L1_cell_probs_data = CSV.read((@__DIR__)*"/../data/L1_cell_probs.csv",DataFrame)
L1_cell_probs_data_lwr = CSV.read((@__DIR__)*"/../data/L1_cell_probs_ci_lwr.csv",DataFrame)
L1_cell_probs_data_upr = CSV.read((@__DIR__)*"/../data/L1_cell_probsci_upr.csv",DataFrame)
lwr = [(L1_cell_probs_data[i,j]) .- (L1_cell_probs_data_lwr[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
upr = [(L1_cell_probs_data_upr[i,j]) .- (L1_cell_probs_data[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:false,:false,:false,:false]
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
colours = [1;3;2;4]
for (c,col) in enumerate(names(L1_cell_probs_data))
    plot!(log10.(1:2:21),L1_cell_probs_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="Cell-wise error")
plot!(legend=:outerbottomright)
plot!()
@add_lines!(L1_cell_probs_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/L1_cell_probs_error_formatted.pdf")
end 
p()

function p()
l2_cdf_data = CSV.read((@__DIR__)*"/../data/L2.csv",DataFrame)
l2_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/L2_ci_lwr.csv",DataFrame)
l2_cdf_data_upr = CSV.read((@__DIR__)*"/../data/L2ci_upr.csv",DataFrame)
markerstyles_vec = [:false,:false,:false,:false]
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
colours = [1;3;2;4]
lwr = [(l2_cdf_data[i,j]) .- (l2_cdf_data_lwr[i,j]) for i in 1:size(l2_cdf_data,1), j in 1:size(l2_cdf_data,2)]
 upr = [(l2_cdf_data_upr[i,j]) .- (l2_cdf_data[i,j]) for i in 1:size(l2_cdf_data,1), j in 1:size(l2_cdf_data,2)]
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
for (c,col) in enumerate(names(l2_cdf_data))
    plot!(log10.(1:2:21),l2_cdf_data[:,col],ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L² error - CDF")
plot!(legend=:outerbottomright)
plot!()
@add_lines!(l2_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l2_cdf_error_formatted.pdf")
end 
p()

function p()
_names = ["DG"; "DG (limit)"; "Unif"; "QBD-RAP"]
l1_cdf_data = CSV.read((@__DIR__)*"/../data/L1.csv",DataFrame)
l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/L1_ci_lwr.csv",DataFrame)
l1_cdf_data_upr = CSV.read((@__DIR__)*"/../data/L1ci_upr.csv",DataFrame)
markerstyles_vec = [:false,:false,:false,:false]
lwr = [(l1_cdf_data[i,j]) .- (l1_cdf_data_lwr[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
 upr = [(l1_cdf_data_upr[i,j]) .- (l1_cdf_data[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
 colours = [1;3;2;4]

@halfwidth_plot_defaults()
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
for (c,col) in enumerate(names(l1_cdf_data))
    plot!(log10.(1:2:21),l1_cdf_data[:,col],ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
plot!(legend=:outerbottomright)
plot!()
@add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end
p()