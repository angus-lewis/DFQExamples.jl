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
function p()
ks_data = CSV.read((@__DIR__)*"/../data/errors/ks_error.csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../data/errors/ks_error_lwr.csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../data/errors/ks_error_upr.csv",DataFrame)

plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
colours = [1;3;2;4]
for (c,col) in enumerate(names(ks_data))
    if c!=2
        plot!(log10.(1:2:21),log10.(ks_data[:,col]), 
            ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error")
plot!(legend=(0.7,0.92))
plot!()
 @add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/ks_error_formatted.pdf")
end
p()

function p()
l1_cdf_data = CSV.read((@__DIR__)*"/../data/errors/l1_error.csv",DataFrame)
l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/errors/l1_error_lwr.csv",DataFrame)
l1_cdf_data_upr = CSV.read((@__DIR__)*"/../data/errors/l1_error_upr.csv",DataFrame)
colours = [1;3;2;4]
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cdf_data))
    if c!=2
        plot!(log10.(1:2:21),log10.(l1_cdf_data[:,col]),
            ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
plot!(legend=(0.7,0.92))
plot!(ylims=(-2.5,-0.55))
plot!()
 @add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end
p()


function p()
l1_cell_probs = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error.csv",DataFrame)
l1_cell_probs_lwr = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error_lwr.csv",DataFrame)
l1_cell_probs_upr = CSV.read((@__DIR__)*"/../data/errors/l1_cell_probs_error_upr.csv",DataFrame)
@halfwidth_plot_defaults()
colours = [1;3;2;4]
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cell_probs))
    if c!=2
        plot!(log10.(1:2:21),log10.(l1_cell_probs[:,col]),
            ribbon=(log10.(l1_cell_probs[:,col])-log10.(l1_cell_probs_lwr[:,col]),log10.(l1_cell_probs_upr[:,col])-log10.(l1_cell_probs[:,col])),
            label=col,linestyle=linestyles_vec[c],
            # marker=markerstyles_vec[c],
            linewidth=2,color=colours[c])
    end
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - cell probabilities")
plot!(legend=(0.7,0.92))
plot!()
 @add_lines!(l1_cell_probs,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cell_probs_error_formatted.pdf")
end
p()

function p()
@fullwidth_plot_defaults()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
dg_data = CSV.read((@__DIR__)*"/../data/first_return_cdf_approximations/order_21_model_dg1.csv",DataFrame)
qbdrap_data = CSV.read((@__DIR__)*"/../data/first_return_cdf_approximations/order_21_model_qbdrap4.csv",DataFrame)
sim_cdf = CSV.read("fluidfluid_discontinuous/data/sims_cdf_evaluated.csv",DataFrame)
plot(dg_data.x,dg_data.phase_4, label="DG",color=1;linewidth=2,linestyle=linestyles_vec[1])
plot!(dg_data.x,qbdrap_data.phase_4, label="QBD-RAP",color=4;linewidth=2,linestyle=linestyles_vec[4])
plot!(dg_data.x,sim_cdf.phase_4, label="Sim",color=5,linewidth=2,linestyle=:dot)
plot!(title="CDF of X at first return in phase 00")
plot!(xlabel="x")
plot!(ylabel="probability")
plot!(legend=(0.15,0.9))
savefig((@__DIR__)*"/phase_4_cdf.pdf")
end
p()
