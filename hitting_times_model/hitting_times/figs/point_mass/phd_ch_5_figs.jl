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
ks_data = CSV.read((@__DIR__)*"/../../data/errors/ks_"*ic*".csv",DataFrame)
ks_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/ks_lwr_"*ic*".csv",DataFrame)
ks_data_upr = CSV.read((@__DIR__)*"/../../data/errors/ks_upr_"*ic*".csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
colours = [1;2;3;4]
for (c,col) in enumerate(names(ks_data))
    plot!(1:2:21,log10.(ks_data[:,col]), 
        ribbon=(log10.(ks_data[:,col])-log10.(ks_data_lwr[:,col]),log10.(ks_data_upr[:,col])-log10.(ks_data[:,col])),
        label=col,linestyle=linestyles_vec[c],
        # marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Order")
plot!(ylabel="log₁₀ error")
plot!(title="KS error between sim and approx. CDF")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/ks_error_formatted.pdf")



l1_cdf_data = CSV.read((@__DIR__)*"/../../data/errors/l1_"*ic*".csv",DataFrame)
l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../../data/errors/l1_lwr_"*ic*".csv",DataFrame)
l1_cdf_data_upr = CSV.read((@__DIR__)*"/../../data/errors/l1_upr_"*ic*".csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cdf_data))
    plot!(1:2:21,log10.(l1_cdf_data[:,col]),
        ribbon=(log10.(l1_cdf_data[:,col])-log10.(l1_cdf_data_lwr[:,col]),log10.(l1_cdf_data_upr[:,col])-log10.(l1_cdf_data[:,col])),
        label=col,linestyle=linestyles_vec[c],
        # marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Order")
plot!(ylabel="log₁₀ error")
plot!(title="L¹ error between sim and approx CDF")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")

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
sim_data = CSV.read((@__DIR__)*"/../../data/point_mass/sims.csv",DataFrame)
data = CSV.read((@__DIR__)*"/../../data/point_mass/order_21_model_dg1.csv",DataFrame)
T_idx = only(findall(data.t.==1.2))
# sim_cdf = ht_cdf(sim_data).(data.t[1:T_idx],1)
linestyles_vec = [:solid,:dash,:dashdot, :dot]
plot(data.t[1:T_idx],data.phase_1[1:T_idx],label="Order 21, DG",color=colours[1],linestyle=linestyles_vec[1],linewidth=2)
data = CSV.read((@__DIR__)*"/../../data/point_mass/order_21_model_qbdrap4.csv",DataFrame)
plot!(data.t[1:T_idx],data.phase_1[1:T_idx],label="Order 21, QBDRAP",color=colours[4],linestyle=linestyles_vec[4],linewidth=2)
plot!(data.t[1:T_idx],sim_cdf,label="Sim",color=5,linestyle=:dash,linewidth=2)
plot!(xlabel="t")
plot!(ylabel="Probability")
plot!(title="First exit time approximations")
plot!(legend=(0.5,0.9))
savefig((@__DIR__)*"/cdf_order21DG_and_sims.pdf")

