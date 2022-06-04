using CSV
using DataFrames
using Plots 
using JLD2
using DiscretisedFluidQueues
import Distributions
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/error_metrics.jl")

@halfwidth_plot_defaults()
ic = "point_mass"
plotlyjs()
function p()
ks_data = CSV.read((@__DIR__)*"/../errors/ks_point_mass.csv",DataFrame)
ks_data.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/ks_error.csv",DataFrame).DG_limit

ks_data_lwr = CSV.read((@__DIR__)*"/../errors/ks_lwr_point_mass.csv",DataFrame)
ks_data_lwr.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/ks_error_lwr.csv",DataFrame).DG_limit

ks_data_upr = CSV.read((@__DIR__)*"/../errors/ks_upr_point_mass.csv",DataFrame)
ks_data_upr.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/ks_error_upr.csv",DataFrame).DG_limit

lwr = [log10(ks_data[i,j]) .- log10(ks_data_lwr[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
upr = [log10(ks_data_upr[i,j]) .- log10(ks_data[i,j]) for i in 1:size(ks_data,1), j in 1:size(ks_data,2)]
ks_data = log10.(ks_data)
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
markerstyles_vec = [:false,:false,:false,:false,:false]
_names = ["DG"; "DG lim"; "Unif"; "QBD-RAP"; "DG lin lim"]
colours = [1;3;2;4;5]
for (c,col) in enumerate(names(ks_data))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),ks_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error")
plot!(title="KS error")
# plot!(legend=(0.225,0.55))
error_ticks!(plot!())
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/ks_error_formatted.pdf")
end
p()

function p()
@halfwidth_plot_defaults()

L1_cell_probs_data = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_point_mass.csv",DataFrame))
L1_cell_probs_data.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_cell_probs_error.csv",DataFrame).DG_limit)

L1_cell_probs_data_lwr = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_lwr_point_mass.csv",DataFrame))
L1_cell_probs_data_lwr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_cell_probs_error_lwr.csv",DataFrame).DG_limit)

L1_cell_probs_data_upr = log10.(CSV.read((@__DIR__)*"/../errors/cell_probs_upr_point_mass.csv",DataFrame))
L1_cell_probs_data_upr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_cell_probs_error_upr.csv",DataFrame).DG_limit)

lwr = [(L1_cell_probs_data[i,j]) .- (L1_cell_probs_data_lwr[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
upr = [(L1_cell_probs_data_upr[i,j]) .- (L1_cell_probs_data[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
plot()
linestyles_vec = [:solid,:dash,:dashdot,:dot,:dash]
markerstyles_vec = [:false,:false,:false,:false,:false]
_names = ["DG"; "DG lim"; "Unif"; "QBD-RAP"; "DG lin lim"]
colours = [1;3;2;4;5]
for (c,col) in enumerate(names(L1_cell_probs_data))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),L1_cell_probs_data[:,col], ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="Cell-wise error")
# plot!(legend=(0.225,0.55))
plot!()
@add_lines!(L1_cell_probs_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/L1_cell_probs_error_formatted.pdf")
end
p()


function p()
l1_cdf_data = log10.(CSV.read((@__DIR__)*"/../errors/l1_point_mass.csv",DataFrame))
l1_cdf_data.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_error.csv",DataFrame).DG_limit)

l1_cdf_data_lwr = log10.(CSV.read((@__DIR__)*"/../errors/l1_lwr_point_mass.csv",DataFrame))
l1_cdf_data_lwr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_error_lwr.csv",DataFrame).DG_limit)

l1_cdf_data_upr = log10.(CSV.read((@__DIR__)*"/../errors/l1_upr_point_mass.csv",DataFrame))
l1_cdf_data_upr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data2/errors/l1_error_upr.csv",DataFrame).DG_limit)

lwr = [(l1_cdf_data[i,j]) .- (l1_cdf_data_lwr[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
upr = [(l1_cdf_data_upr[i,j]) .- (l1_cdf_data[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
@halfwidth_plot_defaults()
plot()
markerstyles_vec = [:false,:false,:false,:false,:false]
_names = ["DG"; "DG lim"; "Unif"; "QBD-RAP"; "DG lin lim"]
colours = [1;3;2;4;5]
linestyles_vec = [:solid,:dash,:dashdot,:dot,:dash]
for (c,col) in enumerate(names(l1_cdf_data))
    os = collect(1:2:21) 
    (c==5)&&(os.+=1)
    plot!(log10.(os),l1_cdf_data[:,col],ribbon=(lwr[:,c],upr[:,c]),
        label=_names[c],linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
# plot!(legend=(0.225,0.55))
plot!()
@add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end
p()


function p()
    @fullwidth_plot_defaults

    dg_lim_lin = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/point_mass/data2/cdfs/order_21_model_limiter.csv", DataFrame)
    plot(dg_lim_lin.x,dg_lim_lin.phase_1,label="DG lim lin",color=4,linewidth=2)

    qbdrap = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/approx_cdf/order_21model_qbdrap4.csv",DataFrame)
    plot!(qbdrap.x,qbdrap.phase_1,label="QBD-RAP",color=5,linewidth=2)

    sims = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/reflecting_model/transient_distribution/point_mass/datat2/sims_cdf_evaluated.csv",DataFrame)
    plot!(sims.x,sims.phase_1,label="Simulation",color=:grey,linewidth=2)

    plot!(title="CDF",xlabel="x",ylabel="Probability",xlims = (1.5,2.6),ylims = (0.4,0.65))


    savefig((@__DIR__)*"/order_21_CDFs_formatted.pdf")
end
p()
