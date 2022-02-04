using CSV
using DataFrames
using Plots 
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")

@halfwidth_plot_defaults()

plotlyjs()
function p()
ks_data = CSV.read((@__DIR__)*"/../data/meshs_ks_func_count_1.csv",DataFrame)

plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
colours = [1;2;4]
for (c,col) in enumerate(names(ks_data))
    plot!(log10.(1:2:21),ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],color=colours[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error")
plot!(legend=(0.7,0.9))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_ks_error_formatted.pdf")
end
p()

function p()
l2_pdf_data = CSV.read((@__DIR__)*"/../data/meshs_l1_cdf_func_count_1.csv",DataFrame)
colours = [1;2;4]

plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
for (c,col) in enumerate(names(l2_pdf_data))
    plot!(log10.(1:2:21),l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],color=colours[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
plot!(legend=(0.8,0.9))
plot!()
@add_lines!(l2_pdf_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_l1_cdf_error_formatted.pdf")
end
p()

function p()
@fullwidth_plot_defaults()

os = 1:2:7
plot(layout=(1,length(os)))
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
xs = range(eps(),1-eps(),length=200)
for (c,o) in enumerate(os)
    model = BoundedFluidQueue(fill(0.0,1,1),[1.0],1.0)
    
    mesh = DGMesh(0.0:1.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
        linestyle=linestyles_vec[1],color=1,
        # marker=markerstyles_vec[1],
        linewidth=2)

    mesh = DGMesh(0.0:(1.0./o):1.0,1)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "Unif." : false),
        linestyle=linestyles_vec[2],color=2,
        # marker=markerstyles_vec[2],
        linewidth=2)

    mesh = FRAPMesh(0.0:1.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,cdf(d0).(xs,1),subplot=c,label=(c==1 ? "QBD-RAP" : false),
        linestyle=linestyles_vec[3],color=4,
        # marker=markerstyles_vec[3],
        linewidth=2)

    plot!(xs,x->Float64(x>=0.5),subplot=c,label=(c==1 ? "truth" : false),
        linestyle=:dot,
        # marker=markerstyles_vec[3],
        color=:grey,
        linewidth=2)

    plot!(subplot=c, yticks=(c==1 ? (0:0.2:1.2) : false),
        ylims = (-0.1,1.1), xticks = 0.25:0.5:0.75, 
        ylabel=(c==1 ? "probability" : false), grid=false, 
        title = string("Dim. ",o),
        legend=(0.93,0.4))
end
plot!() 
savefig((@__DIR__)*"/cdfs_formatted.pdf")
end
p()


function p()
os = 1:2:7
plot(layout=(1,length(os)))
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
xs = range(eps(),1-eps(),length=200)
for (c,o) in enumerate(os)
    model = BoundedFluidQueue(fill(0.0,1,1),[1.0],1.0)
    
    mesh = DGMesh(0.0:1.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
        linestyle=linestyles_vec[1],color=1,
        # marker=markerstyles_vec[1],
        linewidth=2)

    mesh = DGMesh(0.0:(1.0./o):1.0,1)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "Unif." : false),
        linestyle=linestyles_vec[2],color=2,
        # marker=markerstyles_vec[2],
        linewidth=2)

    mesh = FRAPMesh(0.0:1.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = interior_point_mass(0.5,1,dq)    
    plot!(xs,pdf(d0).(xs,1),subplot=c,label=(c==1 ? "QBD-RAP" : false),
        linestyle=linestyles_vec[3],color=4,
        # marker=markerstyles_vec[3],
        linewidth=2)

    plot!(subplot=c, yticks=(c==1 ? (-2:2:6) : false),
        ylims = (-3,7.2), xticks = 0.25:0.5:0.75, 
        ylabel=(c==1 ? "probability" : false), grid=false, 
        title = string("Dim. ",o),
        legend=(0.1,0.8))
end
plot!() 
savefig((@__DIR__)*"/pdfs_formatted.pdf")
end
p()
