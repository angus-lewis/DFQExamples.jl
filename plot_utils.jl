using Plots, GLM
plotlyjs() 
macro halfwidth_plot_defaults()
    return :(default(
        legendfontsize=16,
        titlefontsize=18,
        xlabelfontsize=18,
        ylabelfontsize=18,
        xtickfontsize=18,
        ytickfontsize=18,
        fg_legend=:false,
    ))
end

macro fullwidth_plot_defaults()
    return :(default(
        legendfontsize=12,
        titlefontsize=16,
        xlabelfontsize=12,
        ylabelfontsize=12,
        xtickfontsize=12,
        ytickfontsize=12,
        size=(600, 300),
        fg_legend=:false,
    ))
end

macro fit_error_lm(data,dir)
    data_name = string(data)
    return quote
        fit_error_lm($data_name,$(esc(data)),$(esc(dir)))
    end
end

macro fit_error_lm_log_first(data,dir)
    data_name = string(data)
    return quote
        fit_error_lm_log_first($data_name,$(esc(data)),$(esc(dir)))
    end
end

function fit_error_lm(data_name,data,dir)
    lines = Dict()
    open(dir*"/"*data_name*"_fit.txt","w") do io 
        for n in names(data)
            m = lm([ones(8) log10.(7:2:21)], data[4:end,n])
            println(io,string("Fit log(KS) = β₀ + β₁log(N), for method ",n))
            println(io,m)
            lines[n] = (x=log10.(7:2:21),y=[ones(8) log10.(7:2:21)]*coef(m),slope=coef(m)[2],p=coeftable(m).cols[4][2])
        end
    end
    return lines
end

function fit_error_lm_log_first(data_name,data,dir)
    data2 = copy(data)
    for n in names(data2)
        transform!(data2, n => (x->log10.(x)) => n)
    end
    lines = fit_error_lm(data_name,data2,dir)
    return lines
end

function error_ticks!(p)
    oldticks = yticks(p)[1]
    if minimum(oldticks[1])<-3
        yticks!(p,oldticks[1],string.("1e",round.(oldticks[1],digits=1)))
    else
        yticks!(p,oldticks[1],[string(10.0.^k+0.00001)[1:5] for k in oldticks[1]])
    end
    # plot!(p,yformatter = :scientific)
    ticks = round.(log10.([1;3;9;21]),digits=2)
    xticks!(p,ticks,["1","3","9","21"])
end

macro add_lines!(data,names_set,dir)
    data_name = string(data)
    return :(add_lines!($data_name,$(esc(data)),$(esc(names_set)),$(esc(dir))))
end

function add_lines!(data_name,data,names_set,dir)
    if "naive"∈names(data)
        data = data[:,Not("naive")]
    end
    if data[end,1] < 0.0
        lines = fit_error_lm(data_name,data,dir)
    else
        lines = fit_error_lm_log_first(data_name,data,dir)
    end
    if "Unif"∈names_set
        names_set = (names_set...,"Order","order13")
    end
    if "DG_limiter"∈names_set
        names_set = (names_set...,"dg2","DG (limit)")
    end
    if "QBDRAP"∈names_set
        names_set = (names_set...,"qbdrap4","QBD-RAP","QBD RAP")
    end
    if "DG"∈names_set
        names_set = (names_set...,"dg1")
    end
    plot!(xlims=(xlims(plot!())[1],log10(21)+0.1))#.+(0,0.05))
    biggest = maximum(abs.([lines[i].y[end] for i in names(data)]))
    for col in names(data)
        distances = (lines[col].y[end].-[lines[i].y[end] for i in names(data[:,Not(col)])])./biggest
        ~, idx = findmin(abs.(distances))
        closest = distances[idx]
        if (lines[col].slope<-0.1)&&(lines[col].p<0.005)
            plot!(
                lines[col].x,lines[col].y,linestyle=:dot,label=false,
                color=:black,linewidth=2,
                # series_annotation=[fill("",length(lines[col].x)-1);string("",round(lines[col].slope,digits=2))],
            )
            if (abs(closest)<0.1)&&(closest<=0.0)
                annotate!((lines[col].x[end]+0.09,lines[col].y[end]-0.025*biggest,string(round(lines[col].slope,digits=2))))
            elseif (abs(closest)<0.1)&&(closest>0.0)
                annotate!((lines[col].x[end]+0.09,lines[col].y[end]+0.025*biggest,string(round(lines[col].slope,digits=2))))
            else
                annotate!((lines[col].x[end]+0.09,lines[col].y[end],string(round(lines[col].slope,digits=2))))
            end
        end
    end
    return lines 
end

