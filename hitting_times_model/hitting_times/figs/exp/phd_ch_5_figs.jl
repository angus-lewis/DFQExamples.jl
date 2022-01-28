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
let
ic = "exp"
ks_data = CSV.read((@__DIR__)*"/../../data/errors/ks_"*ic*".csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/ks_lwr_"*ic*".csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../../data/errors/ks_upr_"*ic*".csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
colours = [1;3;2;4]

for (c,col) in enumerate(names(ks_data))
    plot!(1:2:21,log10.(ks_data[:,col]), 
        ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
        label=col,linestyle=linestyles_vec[c],
        # marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="KS error")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/ks_error_formatted.pdf")
end

let
ic = "exp"
l1_cdf_data = CSV.read((@__DIR__)*"/../../data/errors/l1_"*ic*".csv",DataFrame)
l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/l1_lwr_"*ic*".csv",DataFrame)
l1_cdf_data_upr = CSV.read((@__DIR__)*"/../../data/errors/l1_upr_"*ic*".csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
colours = [1;3;2;4]

for (c,col) in enumerate(names(l1_cdf_data))
    plot!(1:2:21,log10.(l1_cdf_data[:,col]),
        ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
        label=col,linestyle=linestyles_vec[c],
        # marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="L¹ error - CDF")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end
let
colours = [1;3;2;4]
@fullwidth_plot_defaults()
function ht_cdf(sims)
    n_sims=length(sims.t)
    function F(x,i)
        i_idx = sims.φ.==i
        Fxi = sum(sims.t[i_idx].<=x)/n_sims
        return Fxi
    end
    return F
end
sim_data = CSV.read((@__DIR__)*"/../../data/exp/sims.csv",DataFrame)
data = CSV.read((@__DIR__)*"/../../data/exp/order_21_model_dg1.csv",DataFrame)
T_idx = only(findall(data.t.==10.0))
sim_cdf = CSV.read((@__DIR__)*"/../../data/exp/sims_cdf_evaluated.csv",DataFrame).phase_2[1:T_idx]
@show maximum(abs.(data.phase_2-sim_cdf))
@show findall(maximum(abs.(data.phase_2-sim_cdf)).==abs.(data.phase_2-sim_cdf))
linestyles_vec = [:solid,:dash,:dashdot, :dot]
plot(data.t[1:T_idx],data.phase_2[1:T_idx],label="DG - Dim 21",color=colours[1],linestyle=linestyles_vec[1],linewidth=2)
data = CSV.read((@__DIR__)*"/../../data/exp/order_21_model_qbdrap4.csv",DataFrame)
@show maximum(abs.(data.phase_2-sim_cdf))
@show findall(maximum(abs.(data.phase_2-sim_cdf)).==abs.(data.phase_2-sim_cdf))
plot!(data.t[1:T_idx],data.phase_2[1:T_idx],label="QBDRAP - Dim 21",color=colours[4],linestyle=linestyles_vec[4],linewidth=2)
plot!(data.t[1:T_idx],sim_cdf,label="Sim",color=5,linestyle=:dash,linewidth=2)
plot!(xlabel="t")
plot!(ylabel="Probability")
plot!(title="CDF of first hitting time in phase 2")
plot!(legend=(0.7,0.35))
savefig((@__DIR__)*"/cdf_order21DG_and_sims.pdf")
end
let
colours = [1;3;2;4]
@fullwidth_plot_defaults()
data_dg = CSV.read((@__DIR__)*"/../../data/exp/order_21_model_dg1.csv",DataFrame)
data_qbdrap = CSV.read((@__DIR__)*"/../../data/exp/order_21_model_qbdrap4.csv",DataFrame)
linestyles_vec = [:solid,:dash,:dashdot, :dot]
h = data.t[2]-data.t[1]
plot((data_dg.t[1:end-1]+data_dg.t[2:end])/2.0,diff(data_dg.phase_1)./h,label="DG - Phase 1",color=colours[1],linestyle=linestyles_vec[1],linewidth=2)
plot!((data_dg.t[1:end-1]+data_dg.t[2:end])/2.0,diff(data_dg.phase_2)./h,label="DG - Phase 2",color=colours[1],linestyle=linestyles_vec[2],linewidth=2)
plot!((data_qbdrap.t[1:end-1]+data_qbdrap.t[2:end])/2.0,diff(data_qbdrap.phase_1)./h,label="QBD-RAP - Phase 1",color=colours[4],linestyle=linestyles_vec[1],linewidth=2)
plot!((data_qbdrap.t[1:end-1]+data_qbdrap.t[2:end])/2.0,diff(data_qbdrap.phase_2)./h,label="QBD-RAP - Phase 2",color=colours[4],linestyle=linestyles_vec[2],linewidth=2)
plot!(xlabel="t")
plot!(ylabel="Probability")
plot!(title="First hitting time PDFs")
plot!(xlims=(0.0,4/3))
plot!(legend=(0.65,0.9))
savefig((@__DIR__)*"/pdf_order21.pdf")
end