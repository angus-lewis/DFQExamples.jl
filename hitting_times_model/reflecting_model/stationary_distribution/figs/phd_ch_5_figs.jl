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
ks_data = CSV.read((@__DIR__)*"/../data/ks_errors.csv",DataFrame)

plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
colours = [1;2;4]
for (c,col) in enumerate(names(ks_data))
    plot!(log10.(1:2:21),ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
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

function p()
l2_pdf_data = CSV.read((@__DIR__)*"/../data/l2_pdf_errors.csv",DataFrame)
colours = [1;2;4]

plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l2_pdf_data))
    plot!(log10.(1:2:21),l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2,color=colours[c])
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L² error - PDF")
plot!(legend=:outerbottomright)
plot!()
@add_lines!(l2_pdf_data,("Unif","QBDRAP"),@__DIR__)

savefig((@__DIR__)*"/l2_pdf_error_formatted.pdf")
end
p()

function p()
l1_cdf_data = CSV.read((@__DIR__)*"/../data/l1_cdf_errors.csv",DataFrame)
@halfwidth_plot_defaults()
colours = [1;2;4]
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l1_cdf_data))
    plot!(log10.(1:2:21),l1_cdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
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
#             approx_pm = [coeffs[1:N₋(model)];coeffs[end-N₊(model)+1:end]]

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
# plot!(xlabel="Dimension", ylabel="log₁₀ Error", 
#     title="Error between cell masses",
#     legend=:outerbottomright)
# savefig((@__DIR__)*"/L1_cell_probs.pdf")

