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
ic = "point_mass"
ks_data = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnov.csv",DataFrame)
ks_data.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/ks_error.csv",DataFrame).DG_limit

ks_data_lwr = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnov_ci_lwr.csv",DataFrame)
ks_data_lwr.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/ks_error_lwr.csv",DataFrame).DG_limit

ks_data_upr = CSV.read((@__DIR__)*"/../data/kolmogorov_smirnovci_upr.csv",DataFrame)
ks_data_upr.DG_lin_lim = CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/ks_error_upr.csv",DataFrame).DG_limit

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
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error - CDF")
# plot!(legend=(0.8,0.97))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/ks_error_formatted.pdf")
end 
p()

function p()
ic = "point_mass"
L1_cell_probs_data = CSV.read((@__DIR__)*"/../data/L1_cell_probs.csv",DataFrame)
L1_cell_probs_data.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_cell_probs_error.csv",DataFrame).DG_limit)

L1_cell_probs_data_lwr = CSV.read((@__DIR__)*"/../data/L1_cell_probs_ci_lwr.csv",DataFrame)
L1_cell_probs_data_lwr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_cell_probs_error_lwr.csv",DataFrame).DG_limit)

L1_cell_probs_data_upr = CSV.read((@__DIR__)*"/../data/L1_cell_probsci_upr.csv",DataFrame)
L1_cell_probs_data_upr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_cell_probs_error_upr.csv",DataFrame).DG_limit)

lwr = [(L1_cell_probs_data[i,j]) .- (L1_cell_probs_data_lwr[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]
upr = [(L1_cell_probs_data_upr[i,j]) .- (L1_cell_probs_data[i,j]) for i in 1:size(L1_cell_probs_data,1), j in 1:size(L1_cell_probs_data,2)]

p = plot()
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
# plot!(legend=(0.7,0.9))
plot!()
@add_lines!(L1_cell_probs_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/L1_cell_probs_error_formatted.pdf")
end 
p()

# function p()
# l2_cdf_data = CSV.read((@__DIR__)*"/../data/L2.csv",DataFrame)
# l2_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/L2_ci_lwr.csv",DataFrame)
# l2_cdf_data_upr = CSV.read((@__DIR__)*"/../data/L2ci_upr.csv",DataFrame)
#  lwr = [(l2_cdf_data[i,j]) .- (l2_cdf_data_lwr[i,j]) for i in 1:size(l2_cdf_data,1), j in 1:size(l2_cdf_data,2)]
#  upr = [(l2_cdf_data_upr[i,j]) .- (l2_cdf_data[i,j]) for i in 1:size(l2_cdf_data,1), j in 1:size(l2_cdf_data,2)]
# plot()
# _names = ["DG"; "DG lim"; "Unif"; "QBD-RAP"]
# markerstyles_vec = [:false,:false,:false,:false]
# colours = [1;3;2;4]
# linestyles_vec = [:solid,:dash,:dashdot, :dot]
# for (c,col) in enumerate(names(l2_cdf_data))
#     plot!(log10.(1:2:21),l2_cdf_data[:,col],ribbon=(lwr[:,c],upr[:,c]),
#         label=_names[c],linestyle=linestyles_vec[c],
#         marker=markerstyles_vec[c],
#         linewidth=2,color=colours[c])
# end
# plot!(xlabel="Dimension")
# plot!(ylabel="Error"); error_ticks!(plot!())
# plot!(title="L?? error - CDF")
# plot!(legend=:outerbottomright)
# plot!()
# @add_lines!(l2_cdf_data,("Unif","QBDRAP"),@__DIR__)

# savefig((@__DIR__)*"/l2_cdf_error_formatted.pdf")
# end
# p()

function p()
ic = "point_mass"

l1_cdf_data = CSV.read((@__DIR__)*"/../data/L1.csv",DataFrame)
l1_cdf_data.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_error.csv",DataFrame).DG_limit)

l1_cdf_data_lwr = CSV.read((@__DIR__)*"/../data/L1_ci_lwr.csv",DataFrame)
l1_cdf_data_lwr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_error_lwr.csv",DataFrame).DG_limit)

l1_cdf_data_upr = CSV.read((@__DIR__)*"/../data/L1ci_upr.csv",DataFrame)
l1_cdf_data_upr.DG_lin_lim = log10.(CSV.read("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_linearDGwLimiter/reflecting_model/data/"*ic*"/data/errors/l1_error_upr.csv",DataFrame).DG_limit)

lwr = [(l1_cdf_data[i,j]) .- (l1_cdf_data_lwr[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
upr = [(l1_cdf_data_upr[i,j]) .- (l1_cdf_data[i,j]) for i in 1:size(l1_cdf_data,1), j in 1:size(l1_cdf_data,2)]
 markerstyles_vec = [:false,:false,:false,:false,:false]

@halfwidth_plot_defaults()
plot()
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
plot!(title="L?? error - CDF")
# plot!(legend=(0.215,0.54))
# plot!(ylims=ylims(plot!()).+(-0.175,00),xlims=xlims(plot!()).+(-0.0375,0))
@add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l1_cdf_error_formatted.pdf")
end 
p()
# @fullwidth_plot_defaults()
# dir = @__DIR__()
# os = 1:2:7
# plot(layout=(4,length(os)))
# linestyles_vec = [:solid,:solid,:solid,:solid]
# markerstyles_vec = [:cross,:diamond,:circle]
# xs = range(eps(),10-eps(),length=200)
# for (c,o) in enumerate(os)
#     model = BoundedFluidQueue(fill(0.0,2,2),[1.0; 0.0],10.0)
    
#     mesh = DGMesh(0.0:10.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
#         coeffs = f["coeff_matrix"][fun_number,c,1]
#         dt = SFMDistribution(coeffs,dq)   
#         plot!(xs,pdf(dt).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
#             linestyle=linestyles_vec[1],ylims = (-0.1,1.1),  yticks=false,
#             # marker=markerstyles_vec[1],
#             linewidth=2,xticks=false)
#     end
#     plot!(subplot=c, yticks=(c==1 ? (0.1:0.4:1.2) : false),
#         ylabel=(c==1 ? "DG" : false), grid=false, 
#         title = string("Dimension ",o))

#     mesh = DGMesh(0.0:(1.0./o):10.0,1)
#     dq = DiscretisedFluidQueue(model,mesh)
#     jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
#         coeffs = f["coeff_matrix"][fun_number,c,2]
#         dt = SFMDistribution(coeffs,dq) 
#         plot!(xs,pdf(dt).(xs,1),subplot=c+4,label=(c==1 ? "unif." : false),
#             linestyle=linestyles_vec[2],ylims = (-0.1,1.1),  yticks=false,
#             # marker=markerstyles_vec[2],
#             linewidth=2,xticks=false)
#     end
#     (c==1)&&plot!(subplot=c+4, yticks=(c==1 ? (0.1:0.4:1.2) : false), 
#         ylabel=(c==1 ? "unif." : false), grid=false)

#     mesh = DGMesh(0.0:10.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
#         coeffs = f["coeff_matrix"][fun_number,c,3]
#         dt = SFMDistribution(coeffs,dq)   
#         plot!(xs,pdf(dt).(xs,1),subplot=c+8,label=(c==1 ? "DG (limiter)" : false),
#             linestyle=linestyles_vec[3],ylims = (-0.1,1.1),  yticks=false,
#             # marker=markerstyles_vec[3],
#             linewidth=2,xticks=false)
#     end
#     (c==1)&&plot!(subplot=c+8, yticks=(c==1 ? (0.1:0.4:1.2) : false),
#         ylims = (-0.1,1.1), 
#         ylabel=(c==1 ? "DG-limit" : false), grid=false)

#     mesh = FRAPMesh(0.0:10.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
#         coeffs = f["coeff_matrix"][fun_number,c,4]
#         dt = SFMDistribution(coeffs,dq)   
#         plot!(xs,pdf(dt).(xs,1),subplot=c+12,label=(c==1 ? "QBD-RAP" : false),
#             linestyle=linestyles_vec[3],ylims = (-0.1,1.1), yticks=false,
#             # marker=markerstyles_vec[3],
#             linewidth=2,xticks=0:2:10, grid=false)
#     end
#     (c==1)&&plot!(subplot=c+12, yticks=(c==1 ? (0.1:0.4:1.2) : false),
#         ylims = (-0.1,1.1), 
#         ylabel=(c==1 ? "QBD-RAP" : false), grid=false)
#     # plot!(xs,x->2*x*(x<0.5) + Float64(x>=0.5),subplot=c,label=(c==1 ? "truth" : false),
#     #     linestyle=:dot,
#     #     # marker=markerstyles_vec[3],
#     #     color=:grey,
#     #     linewidth=2)

#     # plot!(subplot=c, yticks=(c==1 ? (0.1:0.3:1.2) : false),
#     #     ylims = (-0.1,1.1), xticks = 0.0:3:10, 
#     #     ylabel=(c==1 ? "probability" : false), grid=false, 
#     #     title = string("Dimension ",o),
#     #     legend=:outertopright)
#     plot!(legend=false)
# end
# plot!() 
# savefig((@__DIR__)*"/pdf_formatted.pdf")

# @halfwidth_plot_defaults()
# L1_cell_probs_errors = DataFrame(DG = Float64[], Order_1 = Float64[], DG_limiter = Float64[], QBDRAP = Float64[])
# i=1; t=4
# F4(z) = exp(-1*t)*(exp(z)-1)*(z<=t) + (z>t)*(1-exp(-1*t))
# truth = F4.(1:10)-F4.(0:9)
# truth_pm = zeros(3)
# for (c_o,o) in enumerate(1:2:21)
#     row = zeros(4)
#     for (c_m,m) in enumerate(names(L1_cell_probs_errors))
#         jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
#             model = BoundedFluidQueue(fill(0.0,2,2),[1.0; 0.0],10.0)
    
#             mesh = (c_m<4) ? DGMesh(0.0:10.0,o) : FRAPMesh(0.0:10.0,o)
#             dq = DiscretisedFluidQueue(model,mesh)

#             coeffs = f["coeff_matrix"][fun_number,c_o,c_m]
#             dt = SFMDistribution(coeffs,dq) 

#             x_vals = range(0.5,9.5;length=10)

#             approx = sum(cell_probs(dt).(x_vals,(1:2)'),dims=2)
#             approx_pm = [coeffs[1:N???(model)];coeffs[end-N???(model)+1:end]]

#             (o==3)&&display(approx)
#             row[c_m] = log10(sum(abs.(approx-truth)) + sum(abs.(approx_pm-truth_pm)))
#         end
#     end
#     push!(L1_cell_probs_errors,row)
# end

# file = (@__DIR__)*"/../data/meshs_l1_cell_"*"func_count_$(fun_number)"
# CSV.write(file*".csv",L1_cell_probs_errors)
# q = plot()
# linestyles_vec = [:solid,:dash,:dashdot, :dot]
# markerstyles_vec = [:cross,:diamond,:circle,:dot]
# for (c,reconstruction) in enumerate(names(L1_cell_probs_errors))
#     plot!(log10.(1:2:21),L1_cell_probs_errors[:,reconstruction],label=reconstruction,
#         linestyle=linestyles_vec[c],
#         marker=markerstyles_vec[c],
#         linewidth=2)
# end
# plot!(xlabel="Dimension", ylabel="log?????? Error", 
#     title="Error between cell masses",
#     legend=:outerbottomright)
# savefig((@__DIR__)*"/L1_cell_probs.pdf")

