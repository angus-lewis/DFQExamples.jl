using CSV
using DataFrames
using Plots 
using DiscretisedFluidQueues
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")

@halfwidth_plot_defaults()

plotlyjs()

let
ks_data = CSV.read((@__DIR__)*"/../data/ks_cdf_func_count_6.csv",DataFrame)

plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
for (c,col) in enumerate(names(ks_data))
    plot!(log10.(1:2:21),ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error")
plot!(legend=(0.7,0.9))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/ks_error_formatted.pdf")
end

let
l2_pdf_data = CSV.read((@__DIR__)*"/../data/l2_pdf_func_count_6.csv",DataFrame)

plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
for (c,col) in enumerate(names(l2_pdf_data))
    plot!(log10.(1:2:21),l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L² error - PDF")
plot!(legend=(0.7,0.9))
plot!()
@add_lines!(l2_pdf_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/l2_pdf_error_formatted.pdf")
end

let
@fullwidth_plot_defaults()
os = 1:2:7
plot(layout=(1,length(os)),
size=(700, 400))
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
xs = range(eps(),1-eps(),length=200)
for (c,o) in enumerate(os)
    mesh = FRAPMesh(0.0:1.0,o)
    model = BoundedFluidQueue(fill(0.0,1,1),[1.0],1.0)
    dq = DiscretisedFluidQueue(model,mesh)
    d0 = SFMDistribution((x,i)->-6*x^2+6*x,dq,TrapezoidRule;fun_evals=10_001)    
    plot!(xs,pdf(d0,unnormalised_closing_operator_pdf).(xs,1),subplot=c,label=(c==1 ? "unnormalised" : false),
        linestyle=linestyles_vec[1],
        # marker=markerstyles_vec[1],
        linewidth=2)
    plot!(xs,x->-6*x^2+6*x,subplot=c,label=(c==2 ? "true PDF" : false),
        linestyle=linestyles_vec[2],
        # marker=markerstyles_vec[2],
        linewidth=2)
    plot!(xs,pdf(d0,normalised_closing_operator_pdf).(xs,1),subplot=c,label=(c==3 ? "normalised" : false),
        linestyle=linestyles_vec[3],
        # marker=markerstyles_vec[3],
        linewidth=2)
    plot!(subplot=c, yticks=(c==1 ? (0:0.4:1.6) : false),
        ylims = (-0.1,1.75), xticks = 0.25:0.5:0.75, 
        ylabel=(c==1 ? "density" : false), grid=false, 
        title = string("Dim. ",o),
        legend=(0.825,0.3))
end
plot!() 
savefig((@__DIR__)*"/pdfs_formatted.pdf")
end