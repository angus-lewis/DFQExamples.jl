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

ic = "point_mass"
function p()
    ks_data = CSV.read((@__DIR__)*"/../../data/errors/ks_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/ks_"*ic*".csv",DataFrame)
    ks_data.DG_limit2 = tmp.DG_limit
    
    ks_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/ks_lwr_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/ks_lwr_"*ic*".csv",DataFrame)
    ks_data_lwr.DG_limit2 = tmp.DG_limit
    
    ks_data_upr = CSV.read((@__DIR__)*"/../../data/errors/ks_upr_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/ks_upr_"*ic*".csv",DataFrame)
    ks_data_upr.DG_limit2 = tmp.DG_limit
    
    plot()
    linestyles_vec = [:solid,:dash,:dashdot, :dot, :dash]
    markerstyles_vec = [:cross,:diamond,:circle,:dot,:diamond]
    colours = [1;3;2;4;5]
    _names = ["DG", "DG lim", "Unif","QBD-RAP", "DG lin lim"]
    
    for (c,col) in enumerate(names(ks_data))
        if c<5
            plot!(log10.(1:2:21),log10.(ks_data[:,col]), 
                ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
                label=_names[c],linestyle=linestyles_vec[c],
                # marker=markerstyles_vec[c],
                linewidth=2,color=colours[c])
        else
            plot!(log10.((1:2:21).+1),log10.(ks_data[:,col]),
                ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
                label=_names[c],linestyle=linestyles_vec[c],
                # marker=markerstyles_vec[c],
                linewidth=2,color=colours[c])
        end
    end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error")
# plot!(legend=(0.8,0.8))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/ks_error_formatted.pdf")
end
p()


function p()
    l1_cdf_data = CSV.read((@__DIR__)*"/../../data/errors/l1_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/l1_"*ic*".csv",DataFrame)
    l1_cdf_data.DG_limit2 = tmp.DG_limit
    
    l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/l1_lwr_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/l1_lwr_"*ic*".csv",DataFrame)
    l1_cdf_data_lwr.DG_limit2 = tmp.DG_limit
    
    l1_cdf_data_upr = CSV.read((@__DIR__)*"/../../data/errors/l1_upr_"*ic*".csv",DataFrame)
    tmp = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/errors/l1_upr_"*ic*".csv",DataFrame)
    l1_cdf_data_upr.DG_limit2 = tmp.DG_limit
    _names = ["DG", "DG lim", "Unif","QBDRAP", "DG lin lim"]
    
    plot()
    linestyles_vec = [:solid,:dash,:dashdot, :dot, :dash]
    markerstyles_vec = [:cross,:diamond,:circle,:cross, :diamond]
    colours = [1;3;2;4;5]
    
    for (c,col) in enumerate(names(l1_cdf_data))
        if c < 5 
            plot!(log10.(1:2:21),log10.(l1_cdf_data[:,col]),
                ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
                label=_names[c],linestyle=linestyles_vec[c],
                # marker=markerstyles_vec[c],
                linewidth=2,color=colours[c])
        else
            plot!(log10.((1:2:21).+1),log10.(l1_cdf_data[:,col]),
                ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
                label=_names[c],linestyle=linestyles_vec[c],
                # marker=markerstyles_vec[c],
                linewidth=2,color=colours[c])
        end
    end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
# plot!(legend=(0.95,0.95))
# plot!(legend=(0.215,0.535))
# plot!(ylims=ylims(plot!()).+(-0.13,0.0),xlims=xlims(plot!()).+(-0.03,0.0))
@add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end 
p()


function p()
colours = [1;3;2;4;5]

@fullwidth_plot_defaults()
sim_data = CSV.read((@__DIR__)*"/../../data/point_mass/sims.csv",DataFrame)
data = CSV.read((@__DIR__)*"/../../data/point_mass/order_21_model_dg1.csv",DataFrame)
T_idx = only(findall(data.t.==1.2))
sim_cdf = CSV.read((@__DIR__)*"/../../data/point_mass/sims_cdf_evaluated.csv",DataFrame).phase_1[1:T_idx]
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
plot(data.t[1:T_idx],data.phase_1[1:T_idx],label="DG - 21",color=colours[1],linestyle=linestyles_vec[1],linewidth=2)
data = CSV.read((@__DIR__)*"/../../data/point_mass/order_21_model_qbdrap4.csv",DataFrame)
plot!(data.t[1:T_idx],data.phase_1[1:T_idx],label="QBDRAP - 21",color=colours[4],linestyle=linestyles_vec[4],linewidth=2)
data = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/hitting_times/data/point_mass/order_21_model_dg1.csv",DataFrame)
plot!(data.t[1:T_idx],data.phase_1[1:T_idx],label="DG lin lim - 22",color=colours[5],linestyle=linestyles_vec[5],linewidth=2)
plot!(data.t[1:T_idx],sim_cdf,label="Sim",color=6,linestyle=:dashdot,linewidth=2)
plot!(xlabel="t")
plot!(ylabel="Probability")
plot!(title="First hitting time CDF")
plot!(legend=(0.2,0.9))
savefig((@__DIR__)*"/cdf_order21DG_and_sims.pdf")

end
p()
