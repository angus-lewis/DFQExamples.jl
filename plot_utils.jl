using Plots  
plotlyjs() 
macro halfwidth_plot_defaults()
    return :(default(
        legendfontsize=18,
        titlefontsize=18,
        xlabelfontsize=18,
        ylabelfontsize=18,
        xtickfontsize=18,
        ytickfontsize=18,
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
    ))
end