using CSV
using DataFrames
using Plots 
using JLD2
using DiscretisedFluidQueues
import Distributions
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/error_metrics.jl")
fun_number = 6
@halfwidth_plot_defaults()

plotlyjs()
function p()

ks_data = CSV.read((@__DIR__)*"/../data/meshs_ks_func_count_$(fun_number).csv",DataFrame)

ks_data.DG_lin_lim = CSV.read((@__DIR__)*"/../data2/meshs_ks_func_count_$(fun_number).csv",DataFrame).DG_limiter

plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
markerstyles_vec = [:cross,:diamond,:circle,:dot,:diamond]
for (c,col) in enumerate(names(ks_data))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="KS error - CDF")
# plot!(legend=(0.275,0.625))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_ks_error_formatted.pdf")
end
p()

function p()
l2_pdf_data = CSV.read((@__DIR__)*"/../data/meshs_l2_pdf_func_count_$(fun_number).csv",DataFrame)
l2_pdf_data.DG_lin_lim = CSV.read((@__DIR__)*"/../data2/meshs_l2_pdf_func_count_$(fun_number).csv",DataFrame).DG_limiter

plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
markerstyles_vec = [:cross,:diamond,:circle,:cross,:diamond]
for (c,col) in enumerate(names(l2_pdf_data))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L² error - PDF")
# plot!(legend=(0.25,0.625))
# plot!(ylims=ylims(plot!()).+(-0.3,0))
@add_lines!(l2_pdf_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_l2_pdf_error_formatted.pdf")
end 
p()

function p()
l1_cdf_data = CSV.read((@__DIR__)*"/../data/meshs_l1_cdf_func_count_$(fun_number).csv",DataFrame)
l1_cdf_data.DG_lin_lim = CSV.read((@__DIR__)*"/../data2/meshs_l1_cdf_func_count_$(fun_number).csv",DataFrame).DG_limiter

@halfwidth_plot_defaults()
plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
markerstyles_vec = [:cross,:diamond,:circle,:cross,:diamond]
for (c,col) in enumerate(names(l1_cdf_data))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),l1_cdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="Error"); error_ticks!(plot!())
plot!(title="L¹ error - CDF")
# plot!(legend=(0.225,0.625))
plot!()
@add_lines!(l1_cdf_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_l1_cdf_error_formatted.pdf")
end
p()

function p()
@fullwidth_plot_defaults()
dir = @__DIR__()
os = 1:2:7
plot(layout=(5,length(os)))
linestyles_vec = [:solid,:solid,:solid,:solid]
markerstyles_vec = [:cross,:diamond,:circle]
xs = range(eps(),10-eps(),length=200)
for (c,o) in enumerate(os)
    model = BoundedFluidQueue(fill(0.0,2,2),[1.0; 0.0],10.0)
    
    mesh = DGMesh(0.0:10.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
        coeffs = f["coeff_matrix"][fun_number,c,1]
        dt = SFMDistribution(coeffs,dq)   
        plot!(xs,pdf(dt).(xs,1),subplot=c,label=(c==1 ? "DG" : false),
            linestyle=linestyles_vec[1],ylims = (-0.1,1.1),  yticks=false,
            # marker=markerstyles_vec[1],
            linewidth=2,xticks=false)
    end
    plot!(subplot=c, yticks=(c==1 ? (0.1:0.4:1.2) : false),
        ylabel=(c==1 ? "DG" : false), grid=false, 
        title = string("Dim. ",o))

    mesh = DGMesh(0.0:(1.0./o):10.0,1)
    dq = DiscretisedFluidQueue(model,mesh)
    jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
        coeffs = f["coeff_matrix"][fun_number,c,2]
        dt = SFMDistribution(coeffs,dq) 
        plot!(xs,pdf(dt).(xs,1),subplot=c+4,label=(c==1 ? "Unif." : false),
            linestyle=linestyles_vec[2],ylims = (-0.1,1.1),  yticks=false,
            # marker=markerstyles_vec[2],
            linewidth=2,xticks=false)
    end
    (c==1)&&plot!(subplot=c+4, yticks=(c==1 ? (0.1:0.4:1.2) : false), 
        ylabel=(c==1 ? "Unif." : false), grid=false)

    mesh = DGMesh(0.0:10.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
        coeffs = f["coeff_matrix"][fun_number,c,3]
        dt = SFMDistribution(coeffs,dq)   
        plot!(xs,pdf(dt).(xs,1),subplot=c+8,label=(c==1 ? "DG (limiter)" : false),
            linestyle=linestyles_vec[3],ylims = (-0.1,1.1),  yticks=false,
            # marker=markerstyles_vec[3],
            linewidth=2,xticks=false)
    end
    (c==1)&&plot!(subplot=c+8, yticks=(c==1 ? (0.1:0.4:1.2) : false),
        ylims = (-0.1,1.1), 
        ylabel=(c==1 ? "DG-limit" : false), grid=false)

    mesh = FRAPMesh(0.0:10.0,o)
    dq = DiscretisedFluidQueue(model,mesh)
    jldopen(dir*"/../../../coeffs_matrix.jld2") do f 
        coeffs = f["coeff_matrix"][fun_number,c,4]
        dt = SFMDistribution(coeffs,dq)   
        plot!(xs,pdf(dt).(xs,1),subplot=c+12,label=(c==1 ? "Q-R" : false),
            linestyle=linestyles_vec[3],ylims = (-0.1,1.1), yticks=false,
            # marker=markerstyles_vec[3],
            linewidth=2,xticks=0:2:10, grid=false)
    end
    (c==1)&&plot!(subplot=c+12, yticks=(c==1 ? (0.1:0.4:1.2) : false),
        ylims = (-0.1,1.1), 
        ylabel=(c==1 ? "Q-R" : false), grid=false)
    # plot!(xs,x->2*x*(x<0.5) + Float64(x>=0.5),subplot=c,label=(c==1 ? "truth" : false),
    #     linestyle=:dot,
    #     # marker=markerstyles_vec[3],
    #     color=:grey,
    #     linewidth=2)

    # plot!(subplot=c, yticks=(c==1 ? (0.1:0.3:1.2) : false),
    #     ylims = (-0.1,1.1), xticks = 0.0:3:10, 
    #     ylabel=(c==1 ? "probability" : false), grid=false, 
    #     title = string("Dim. ",o),
    #     legend=:outertopright)
    plot!(legend=false)

    model = BoundedFluidQueue([0 0; 1.0 -1.0],[1.0; 0.0],10.0)

    mesh = DGMesh(0.0:(2/(o+1)):10.0,2)
    dq = DiscretisedFluidQueue(model,mesh)
    thepdf0(z,i) = (norm_pdf(z)/(1-norm_cdf(0)))*(i==1)*(z>=0.0)
    d0 = SFMDistribution(thepdf0,dq)
    dt = integrate_time(d0,build_full_generator(dq),4.0, StableRK4(0.005); limiter=GeneralisedMUSCL)
 
    plot!(xs,pdf(dt).(xs,1),subplot=c+16,label=(c==1 ? "DG lin lim" : false),
        linestyle=linestyles_vec[3],ylims = (-0.1,1.1), yticks=false,
        # marker=markerstyles_vec[3],
        linewidth=2,xticks=0:2:10, grid=false)
    (c==1)&&plot!(subplot=c+16, yticks=(c==1 ? (0.1:0.4:1.2) : false),
        ylims = (-0.1,1.1), 
        ylabel=(c==1 ? "DG lin lim" : false), grid=false)
    plot!(legend=false)
end
plot!() 
savefig((@__DIR__)*"/pdf_formatted.pdf")
end 
p()

function p()
@halfwidth_plot_defaults()
L1_cell_probs_errors = DataFrame(DG = Float64[], Unif = Float64[], DG_limiter = Float64[], QBDRAP = Float64[], DG_lin_lim = Float64[])
# truth = zeros(10)
# truth[5] = 1.0
norm_cdf(x) = Distributions.cdf(Distributions.Normal(2.5,0.5),x)
i=1
F4(z) = (norm_cdf(z-4)-norm_cdf(0))/(norm_cdf(10)-norm_cdf(0))*(i==1)*(z>=4)*(z<10) + (z>=10)
truth = F4.(1:10)-F4.(0:9)
truth_pm = zeros(3)
for (c_o,o) in enumerate(1:2:21)
    row = zeros(5)
    for (c_m,m) in enumerate(names(L1_cell_probs_errors))
        (c_m<5)&&jldopen((@__DIR__)*"/../../../coeffs_matrix.jld2") do f 
            model = BoundedFluidQueue(fill(0.0,2,2),[1.0; 0.0],10.0)
    
            mesh = (c_m<4) ? DGMesh(0.0:10.0,o) : FRAPMesh(0.0:10.0,o)
            dq = DiscretisedFluidQueue(model,mesh)

            coeffs = f["coeff_matrix"][fun_number,c_o,c_m]
            dt = SFMDistribution(coeffs,dq) 

            x_vals = range(0.5,9.5;length=10)

            approx = sum(cell_probs(dt).(x_vals,(1:2)'),dims=2)
            approx_pm = [coeffs[1:N₋(model)];coeffs[end-N₊(model)+1:end]]

            row[c_m] = log10(sum(abs.(approx-truth)) + sum(abs.(approx_pm-truth_pm)))
        end
        (c_m==5)&&begin 
            model = BoundedFluidQueue([0 0; 1.0 -1.0],[1.0; 0.0],10.0)

            mesh = DGMesh(0.0:(2/(o+1)):10.0,2)
            dq = DiscretisedFluidQueue(model,mesh)
            thepdf0(z,i) = (norm_pdf(z)/(1-norm_cdf(0)))*(i==1)*(z>=0.0)
            d0 = SFMDistribution(thepdf0,dq)
            dt = integrate_time(d0,build_full_generator(dq),4.0, StableRK4(0.005); limiter=GeneralisedMUSCL)
            coeffs = dt.coeffs

            x_vals = (mesh.nodes[1:end-1]+mesh.nodes[2:end])./2

            approx = sum(reshape(cell_probs(dt).(x_vals,1),Int((o+1)/2),10),dims=1)[:]
            approx_pm = [coeffs[1:N₋(model)];coeffs[end-N₊(model)+1:end]]

            row[c_m] = log10(sum(abs.(approx-truth)) + sum(abs.(approx_pm-truth_pm)))
        end 
    end
    push!(L1_cell_probs_errors,row)
end

file = (@__DIR__)*"/../data/meshs_l1_cell_"*"func_count_$(fun_number)"
CSV.write(file*".csv",L1_cell_probs_errors)
q = plot()
linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
markerstyles_vec = [:cross,:diamond,:circle,:dot,:diamond]
for (c,reconstruction) in enumerate(names(L1_cell_probs_errors))
    os = collect(1:2:21)
    (c==5)&&(os.+=1)
    plot!(log10.(os),L1_cell_probs_errors[:,reconstruction],label=reconstruction,
        linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension", ylabel="Error", 
    title="Error between cell probabilities",)
    # legend=(0.225,0.625))
error_ticks!(plot!())
@add_lines!(L1_cell_probs_errors,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/L1_cell_probs.pdf")
end
p()

function p() 
    @fullwidth_plot_defaults 
    model = BoundedFluidQueue([0 0; 1.0 -1.0],[1.0; 0.0],10.0)

    mesh = DGMesh(0.0:(2/(21+1)):10.0,2)
    dq = DiscretisedFluidQueue(model,mesh)
    thepdf0(z,i) = (norm_pdf(z)/(1-norm_cdf(0)))*(i==1)*(z>=0.0)
    d0 = SFMDistribution(thepdf0,dq)
    dt = integrate_time(d0,build_full_generator(dq),4.0, StableRK4(0.005); limiter=GeneralisedMUSCL)

    x_vals = range(5,8,length=201)

    dg_lin_lim = pdf(dt).(x_vals,1)

    frapmesh = FRAPMesh(0.0:1.0:10.0,21)
    frapdq = DiscretisedFluidQueue(model,frapmesh)
    frapd0 = SFMDistribution(thepdf0,frapdq)
    frapdt = integrate_time(frapd0,build_full_generator(frapdq),4.0, StableRK4(0.005))

    frap = pdf(frapdt).(x_vals,1)

    # linestyles_vec = [:solid,:dash,:dashdot, :dot,:dash]
    # markerstyles_vec = [:cross,:diamond,:circle,:dot,:diamond]

    plot(x_vals,x->thepdf0(x-4.0,1),label="true PDF",linestyle=:solid,color=:grey)
    plot!(x_vals,frap,label="QBDRAP",linestyle=:solid,color=4,linewidth=2)
    plot!(x_vals,dg_lin_lim,label="DG-lin-lim",linestyle=:solid,color=5,linewidth=2)
    plot!(xlabel="x",ylabel="Density",title="PDF")
    savefig((@__DIR__)*"/order_21_PDFs_formatted.pdf")
end 
p()