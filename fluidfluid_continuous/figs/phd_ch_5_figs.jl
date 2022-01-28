using CSV
using DataFrames
using Plots 
# using PlotlyJS
using JLD2
using DiscretisedFluidQueues
import Distributions
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/error_metrics.jl")

@halfwidth_plot_defaults()

plotlyjs()

ks_data = CSV.read((@__DIR__)*"/../data/errors/ks_error.csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../data/errors/ks_error_lwr.csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../data/errors/ks_error_upr.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
colours = [1;3;2;4]
for (c,col) in enumerate(names(ks_data))
    if c!=2
        plot!(1:2:21,log10.(ks_data[:,col]), 
            ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="KS error")
plot!(legend=(0.7,0.9))
plot!()
savefig((@__DIR__)*"/ks_error_formatted.pdf")



l1_cdf_data = CSV.read((@__DIR__)*"/../data/errors/l1_error.csv",DataFrame)
l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/errors/l1_error_lwr.csv",DataFrame)
l1_cdf_data_upr = CSV.read((@__DIR__)*"/../data/errors/l1_error_upr.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cdf_data))
    if c!=2
        plot!(1:2:21,log10.(l1_cdf_data[:,col]),
            ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="L¹ error - CDF")
plot!(legend=(0.7,0.9))
plot!()
savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")


l1_cell_probs = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error.csv",DataFrame)
l1_cell_probs_lwr = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error_lwr.csv",DataFrame)
l1_cell_probs_upr = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error_upr.csv",DataFrame)
@halfwidth_plot_defaults()
p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cell_probs))
    if c!=2
        plot!(1:2:21,log10.(l1_cell_probs[:,col]),
            ribbon=(log10.(l1_cell_probs[:,col])-log10.(l1_cell_probs_lwr[:,col]),log10.(l1_cell_probs_upr[:,col])-log10.(l1_cell_probs[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="L¹ error - cell probabilities")
plot!(legend=(0.7,0.9))
plot!()
savefig((@__DIR__)*"/l1_cell_probs_error_formatted.pdf")
