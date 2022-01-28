using CSV
using DataFrames
using Plots 
using DiscretisedFluidQueues
include("/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/plot_utils.jl")

@halfwidth_plot_defaults()

plotlyjs()

ks_data = CSV.read((@__DIR__)*"/../data/ks_cdf_func_count_7.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
for (c,col) in enumerate(names(ks_data))
    plot!(1:2:21,ks_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="KS error")
plot!(legend=(0.7,0.9))
plot!()
savefig((@__DIR__)*"/ks_error_formatted.pdf")



l2_pdf_data = CSV.read((@__DIR__)*"/../data/l2_pdf_func_count_7.csv",DataFrame)

p = plot()
linestyles_vec = [:solid,:dash,:dashdot]
markerstyles_vec = [:cross,:diamond,:circle]
for (c,col) in enumerate(names(l2_pdf_data))
    plot!(1:2:21,l2_pdf_data[:,col],
        label=col,linestyle=linestyles_vec[c],
        marker=markerstyles_vec[c],
        linewidth=2)
end
plot!(xlabel="Dimension")
plot!(ylabel="log₁₀ error")
plot!(title="L² error - PDF")
plot!(legend=(0.7,0.9))
plot!()
savefig((@__DIR__)*"/l2_pdf_error_formatted.pdf")

