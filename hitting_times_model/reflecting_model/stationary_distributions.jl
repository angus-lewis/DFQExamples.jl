include((@__DIR__)*"/model_def.jl")
orders = 1:2:5

for (c_o,o) in enumerate(orders)
    for (c_m,m) in enumerate((:dg,:order1,:qbdrap))
        @eval $(Symbol(m,"_d0",o)),$(Symbol(m,"_B",o)),$(Symbol(m,"_dq",o)) = $(Symbol(m,"_model"))($o)
        @eval $(Symbol(m,"_B",o))[:,3] .= 1.0
        @eval bvec = zeros(1,size($(Symbol(m,"_B",o)),1))
        bvec[3] = 1.0
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
    -eps(), 10.0+eps();
    # ylim=(-0.1,2.1),
    xticks=false,
    yticks=false,
    grid=false,
    label=false,
    kwargs...,
)

true_pdf, true_cdf = stationary_distribution_x(model)

for (c_o,o) in enumerate(orders)
    for (c_m,m) in enumerate((:dg,:order1,:qbdrap))
        for i in 1:2
            @eval std_plot!($(Symbol("p",1)).layout.grid[$c_m,$c_o],x->$(Symbol(m,"_pdf",o))(x,$i))
            @eval std_plot!($(Symbol("p",1)).layout.grid[$c_m,$c_o],x->true_pdf(x)[$i])
            @eval std_plot!($(Symbol("q",1)).layout.grid[$c_m,$c_o],x->$(Symbol(m,"_cdf",o))(x,$i))
            @eval std_plot!($(Symbol("q",1)).layout.grid[$c_m,$c_o],x->true_cdf(x)[$i])
        end
        for plt in (:p,)#,:q) 
            if c_o==1
                if plt==:p
                    @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],yticks=0:0.02:0.06)
                else 
                    @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],yticks=0:0.01:0.05)
                end
            end
            if c_m==length(orders)
                @eval plot!($(Symbol(plt,1)).layout.grid[$c_m,$c_o],xticks=0:2:10)
            end
        end
    end
end
plot!(p1)
plot!(q1)
