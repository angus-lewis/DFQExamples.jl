using CSV
using DataFrames
using Plots 
using JLD2
using DiscretisedFluidQueues
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")

@halfwidth_plot_defaults()

plotlyjs()

ks_data = CSV.read((@__DIR__)*"/../data/meshs_ks_func_count_4.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:dot]
for (c,col) in enumerate(names(ks_data))
    plot!(1:2:21,ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Order")
plot!(ylabel="log₁₀ error")
plot!(title="KS error between true CDF and reconstruction")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/meshs_ks_error_formatted.pdf")



l2_pdf_data = CSV.read((@__DIR__)*"/../data/meshs_l2_pdf_func_count_4.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot]
markerstyles_vec = [:cross,:diamond,:circle,:cross]
for (c,col) in enumerate(names(l2_pdf_data))
    plot!(1:2:21,l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Order")
plot!(ylabel="log₁₀ error")
plot!(title="L² error between true CDF and reconstruction")
plot!(legend=:outerbottomright)
plot!()
savefig((@__DIR__)*"/meshs_l2_pdf_error_formatted.pdf")

@fullwidth_plot_defaults()
dir = @__DIR__()
jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
    f["coeff_matrix"][1,1,1].coeffs
end
# os = 1:2:7
# plot(layout=(1,length(os)))
# linestyles_vec = [:solid,:dash,:dashdot]
# markerstyles_vec = [:cross,:diamond,:circle]
# xs = range(eps(),1-eps(),length=200)
# for (c,o) in enumerate(os)
#     model = BoundedFluidQueue(fill(0.0,1,1),[1.0],1.0)
    
#     mesh = DGMesh(0.0:1.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)      
#     plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
#         linestyle=linestyles_vec[1],
#         # marker=markerstyles_vec[1],
#         linewidth=2)

#     mesh = DGMesh(0.0:(1.0./o):1.0,1)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)    
#     plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "unif." : false),
#         linestyle=linestyles_vec[2],
#         # marker=markerstyles_vec[2],
#         linewidth=2)

#     mesh = FRAPMesh(0.0:1.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)    
#     plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "QBD-RAP" : false),
#         linestyle=linestyles_vec[3],
#         # marker=markerstyles_vec[3],
#         linewidth=2)

#     plot!(xs,x->2*x*(x<0.5) + Float64(x>=0.5),subplot=c,label=(c==1 ? "truth" : false),
#         linestyle=:dot,
#         # marker=markerstyles_vec[3],
#         color=:grey,
#         linewidth=2)

#     plot!(subplot=c, yticks=(c==1 ? (0:0.2:1.2) : false),
#         ylims = (-0.1,1.1), xticks = 0.25:0.5:0.75, 
#         ylabel=(c==1 ? "probability" : false), grid=false, 
#         title = string("Order ",o),
#         legend=(0.93,0.4))
# end
# plot!() 
# savefig((@__DIR__)*"/cdfs_formatted.pdf")




# os = 1:2:7
# plot(layout=(1,length(os)))
# linestyles_vec = [:solid,:dash,:dashdot]
# markerstyles_vec = [:cross,:diamond,:circle]
# xs = range(eps(),1-eps(),length=200)
# for (c,o) in enumerate(os)
#     model = BoundedFluidQueue(fill(0.0,1,1),[1.0],1.0)
    
#     mesh = DGMesh(0.0:1.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)    
#     plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
#         linestyle=linestyles_vec[1],
#         # marker=markerstyles_vec[1],
#         linewidth=2)

#     mesh = DGMesh(0.0:(1.0./o):1.0,1)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)    
#     plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "unif." : false),
#         linestyle=linestyles_vec[2],
#         # marker=markerstyles_vec[2],
#         linewidth=2)

#     mesh = FRAPMesh(0.0:1.0,o)
#     dq = DiscretisedFluidQueue(model,mesh)
#     d0 = SFMDistribution((x,i)->2*Float64(x<=0.5),dq,TrapezoidRule;fun_evals=10_001)    
#     plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "QBD-RAP" : false),
#         linestyle=linestyles_vec[3],
#         # marker=markerstyles_vec[3],
#         linewidth=2)

#     plot!(xs,x->2*Float64(x<=0.5),subplot=c,label=(c==1 ? "truth" : false),
#         linestyle=:dot,
#         # marker=markerstyles_vec[3],
#         color=:grey,
#         linewidth=2)

#     plot!(subplot=c, yticks=(c==1 ? (-2:2:6) : false),
#         ylims = (-3,7.2), xticks = 0.25:0.5:0.75, 
#         ylabel=(c==1 ? "probability" : false), grid=false, 
#         title = string("Order ",o),
#         legend=(0.1,0.8))
# end
# plot!() 
# savefig((@__DIR__)*"/pdfs_formatted.pdf")
