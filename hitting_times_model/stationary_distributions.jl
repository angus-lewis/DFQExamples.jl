## Run this from the hitting_times_model directory
model_str = (@__DIR__)*"/absorbing_model/"
pth = model_str*"stationary_distribution/"
include(model_str*"model_def.jl")
mkpath(pth*"data")
mkpath(pth*"figs")

include((@__DIR__)*"default_params.jl")
# playing with meta programming, sorry its hard to read...
for (c_o,o) in enumerate(orders)
    for (c_m,m) in enumerate((:dg,:order1,:qbdrap))
        @eval $(Symbol(m,"_d0",o)),$(Symbol(m,"_B",o)),$(Symbol(m,"_dq",o)) = $(Symbol(m,"_model"))($o)
        @eval $(Symbol(m,"_B",o))[:,1] .= 1.0
        @eval bvec = zeros(1,size($(Symbol(m,"_B",o)),1))
        bvec[1] = 1.0
        @eval $(Symbol(m,"_π",o)) = SFMDistribution(bvec/Matrix($(Symbol(m,"_B",o))),$(Symbol(m,"_dq",o)))
        @eval $(Symbol(m,"_pdf",o)) = pdf($(Symbol(m,"_π",o)))
        @eval $(Symbol(m,"_cdf",o)) = cdf($(Symbol(m,"_π",o)))
    end
end

p1 = plot(layout = (3,length(orders))) 
q1 = plot(layout = (3,length(orders)))

std_plot!(p,fun,args...; kwargs...) = plot!(
    p,
    fun,
    -eps(), model.b+eps();
    # ylim=(-0.1,2.1),
    xticks=false,
    yticks=false,
    grid=false,
    label=false,
    kwargs...,
)

_pdf, _cdf = sfm.stationary_distribution_x(model) 
true_pdf(x,i) = _pdf(x)[i]
true_cdf(x,i) = _cdf(x)[i]

L1_cdf_errors = DataFrame(DG = Float64[], Order_1 = Float64[], QBDRAP = Float64[])
ks_errors = DataFrame(DG = Float64[], Order_1 = Float64[], QBDRAP = Float64[])
L1_pdf_errors = DataFrame(DG = Float64[], Order_1 = Float64[], QBDRAP = Float64[])
L2_pdf_errors = DataFrame(DG = Float64[], Order_1 = Float64[], QBDRAP = Float64[])

L1_cdf_error_row = zeros(1,3)
ks_error_row = zeros(1,3)
L1_pdf_error_row = zeros(1,3)
L2_pdf_error_row = zeros(1,3)

# what a truely terrible piece of code! 
for (c_o,o) in enumerate(orders)
    for (c_m,m) in enumerate((:dg,:order1,:qbdrap))
        for i in 1:2
            # @eval pm = sum($(Symbol(m,"_π",o))[[1:2;end-1:end]])
            @eval std_plot!($(Symbol("p",1)).layout.grid[$c_m,$c_o],x->$(Symbol(m,"_pdf",o))(x,$i))
            @eval std_plot!($(Symbol("p",1)).layout.grid[$c_m,$c_o],x->true_pdf(x,$i))
            @eval std_plot!($(Symbol("q",1)).layout.grid[$c_m,$c_o],x->$(Symbol(m,"_cdf",o))(x,$i))
            @eval std_plot!($(Symbol("q",1)).layout.grid[$c_m,$c_o],x->true_cdf(x,$i))

        end
        for plt in (:p,:q) 
            if c_o==1
                @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],ylabel=$(uppercase(string(m))))
                if plt==:p
                    @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],yticks=0:0.002:0.006)
                else 
                    @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],yticks=0:0.01:0.05)
                end
            end
            if c_m==1
                @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],title=$(string(o)))
            end
            if c_m==3
                @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],xticks=0:5:10)
            end
        end
        @eval (L1_cdf_error_row[$c_m] = log10(DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_cdf",o))(x,1),x->true_cdf(x,1),range(-eps(),10+eps(),length=n_err_evals)) + 
                                                        DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_cdf",o))(x,2),x->true_cdf(x,2),range(-eps(),10+eps(),length=n_err_evals))))
        @eval (ks_error_row[$c_m] = log10(DiscretisedFluidQueues.kolmogorov_smirnov(x->$(Symbol(m,"_cdf",o))(x,1),x->true_cdf(x,1),range(-eps(),10+eps(),length=n_err_evals)) + 
                                                        DiscretisedFluidQueues.kolmogorov_smirnov(x->$(Symbol(m,"_cdf",o))(x,2),x->true_cdf(x,2),range(-eps(),10+eps(),length=n_err_evals))))
        @eval (L1_pdf_error_row[$c_m] = log10(DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_pdf",o))(x,1),x->true_pdf(x,1),range(-eps(),10+eps(),length=n_err_evals)) + 
                                                        DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_pdf",o))(x,2),x->true_pdf(x,2),range(-eps(),10+eps(),length=n_err_evals))))
        @eval (L2_pdf_error_row[$c_m] = log10(DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_pdf",o))(x,1),x->true_pdf(x,1),range(-eps(),10+eps(),length=n_err_evals),2) + 
                                                        DiscretisedFluidQueues.Lp(x->$(Symbol(m,"_pdf",o))(x,2),x->true_pdf(x,2),range(-eps(),10+eps(),length=n_err_evals),2)))
    end
    push!(L1_cdf_errors,L1_cdf_error_row)
    push!(ks_errors,ks_error_row)
    push!(L1_pdf_errors,L1_pdf_error_row)
    push!(L2_pdf_errors,L2_pdf_error_row)
end
file = pth*"figs/"
savefig(p1,file*"pdfs.svg")
savefig(q1,file*"cdfs.svg")
titles = [
    "L¹ error between the true CDF and approximations";
    "KS error between the true CDF and approximations";
    "L¹ error between the true PDF and approximations";
    "L² error between the true PDF and approximations";
]
for (c_err,err) in enumerate((:L1_cdf_errors,:ks_errors,:L1_pdf_errors,:L2_pdf_errors))
    file = pth*"data/"*string(err)
    @eval CSV.write($file*".csv",$err)
    @eval q = plot()
    df = @eval $err
    for m in names(df)
        @eval plot!(q,orders,$(err)[:,$(string(m))],label=$m)
    end
    @eval plot!(q,xlabel="Order", ylabel="log₁₀ Error", 
        title=titles[$c_err],
        legend=:outertopright)
    #display(q)
    # #display((@__DIR__)*)
    file = pth*"figs/"*string(err)
    @eval savefig(q,$file*".svg")
end 
