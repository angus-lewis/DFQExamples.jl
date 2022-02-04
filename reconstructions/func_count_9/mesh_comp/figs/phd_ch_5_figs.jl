using CSV
using DataFrames
using Plots 
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")

@halfwidth_plot_defaults()

plotlyjs()
function p()
ks_data = CSV.read((@__DIR__)*"/../data/meshs_ks_func_count_9.csv",DataFrame)

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
plot!(legend=(0.2,0.5))
plot!()
@add_lines!(ks_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_ks_error_formatted.pdf")
end
p()

function p()
l2_pdf_data = CSV.read((@__DIR__)*"/../data/meshs_l2_pdf_func_count_9.csv",DataFrame)

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
plot!(title="L² error - PDF")
plot!(legend=(0.15,0.5))
plot!()
@add_lines!(l2_pdf_data,("Unif","QBDRAP"),@__DIR__)
savefig((@__DIR__)*"/meshs_l2_pdf_error_formatted.pdf")
end
p()
