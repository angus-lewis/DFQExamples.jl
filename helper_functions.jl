function make_approximations(orders,approx_types,d0_map,a_map=identity)
    models = Dict{Any,Dict{Any,Dict}}()
    for (c_o,o) in enumerate(orders)
        tmp_dict = Dict()
        dg_counter = 0
        for (c_m,m) in enumerate(approx_types)

            ~,B,dq = @eval $(Symbol(m,"_model"))($o)

            d0 = d0_map(dq)

            (string(m)=="dg")&&(dg_counter+=1)
            (dg_counter==2) ? (l=true) : (l=false)

            dt = a_map(d0,B,m,transient_time,0.005,l)

            the_pdf = pdf(dt)
            the_cdf = cdf(dt)
            tmp_dict[string(m,c_m)] = Dict(
                "pdf"=>the_pdf,
                "cdf"=>the_cdf,
                "coeffs_mapped"=>dt,
                "coeffs_0"=>d0,
                "dq"=>dq,
                "B"=>B,
            )
        end
        models[o] = tmp_dict
    end
    return models
end 

macro evaluate!(models,funs,grid,_phases)
    expr = quote 
        evaluate!(
            $(esc(models)),
            $(esc(funs)),
            $(esc(grid)),
            $(esc(_phases)),
            $(string.(esc(funs).args[1].args)),
        )
    end
    return expr
end

function evaluate!(models,funs::Tuple{Functions},grid::AbstractVector,
    _phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple})

    orders = sort([l for l in keys(models)])
    k = keys(models[orders[1]])
    model = models[orders[1]][k[1]]["coeffs_mapped"].dq.model
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(sort([i for i in k]))
            for (c_f,f) in enumerate(funs)
                temp_fun = f(models[o][m]["coeffs_mapped"])
                models[c_o][c_m][fun_names[c_f]*"_evaluated"] = temp_fun.(grid[:],_phases)
                models[c_o][c_m]["boundary_masses"] = 
                    models[o][m]["coeffs_mapped"][[1:N₋(model);end-N₊(model)+1:end]]
            end
        end
    end
    return nothing
end
# evaluate!(models,fun::Functions,grid::AbstractVector,_phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple}) = 
#     evaluate!(models,(fun,),grid,_phases,fun_names)

macro evaluate!(models,funs,grid)
    expr = quote 
        evaluate!(
            $(esc(models)),
            $(esc(funs)),
            $(esc(grid)),
            $(string.(esc(funs).args[1].args)),
        )
    end
    return expr
end
function evaluate!(models,funs::Tuple{Functions},grid::AbstractVector,fun_names::Union{AbstractVector{String},Tuple}) 
    o = sort([l for l in keys(models)])
    k = keys(models[o[1]])
    model = models[o[1]][k[1]]["coeffs_mapped"].dq.model
    evaluate!(models,funs,grid,1:n_phases(model),fun_names)
    return nothing
end
# evaluate!(models,fun::Function,grid::AbstractVector,fun_names::Union{AbstractVector{String},Tuple}) = 
#     evaluate!(models,(fun,),grid,fun_names)

macro evaluate(sims,funs,grid,_phases)
    return quote 
        evaluate($(esc(sims)),$(esc(funs)),$(esc(grid)),$(esc(_phases)),$(string.(esc(funs).args[1].args)))
    end
end
function evaluate(sims::Simulation,funs::Tuple{Functions},grid::AbstractVector,_phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple})
    evaluated = Dict{String,AbstractArray}()
    for (c_f,f) in enumerate(funs)
        temp_fun = f(sims)
        evaluated[fun_names[c_f]*"_evaluated"] = temp_fun.(grid,_phases)
        evaluated["boundary_masses"] = sim_point_masses(sims)
    end
    return evaluated
end
macro evaluate(sims,funs,grid)
    return quote 
        evaluate($(esc(sims)),$(esc(funs)),$(esc(grid)),$(string.(esc(funs).args[1].args)))
    end
end
function evaluate(sims::Simulation,funs::Tuple{Functions},grid::AbstractVector,fun_names::Union{AbstractVector{String},Tuple})
    return evaluate(sims,funs,grid,1:n_phases(sims.model),fun_names)
end

function simulate_model(model,stopping_time,n_sims,init,pth,rng::Random.AbstractRNG=Randon.default_rng())
    SFM0 = init(n_sims,rng)
    sims = simulate(model,stopping_time,SFM0,rng)
    df = DataFrame(t=sims.t,φ=sims.φ,X=sims.X,n=sims.n)
    CSV.write(pth*"/data/sim.csv",df)
    return sims
end

function compute_errors(models,truth,err_fun,distribution_function_model,distribution_function_sim,grid,_phases)
    orders = sort([l for l in keys(models)])
    k = keys(models[orders[1]])
    err_data = Matrix{Union{Missing, Float64}}(undef,length(orders),length(k))
    evaluated, pm_evaluated = _eval_sim(truth,distribution_function_sim,grid,length(truth.t))
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(sort([i for i in k]))
            err_data[c_o,c_m] = _compute_errors(models[o][m],evaluated,pm_evaluated,err_fun,distribution_function_model,grid,_phases)
        end
    end
    return err_data
end

function _compute_errors(model,evaluated,pm_evaluated,err_fun,distribution_function,grid,_phases)
    return log10(err_fun(distribution_function(model["coeffs_mapped"]),evaluated,pm_evaluated,grid,_phases))
end

function sim_point_masses(sims)
    pm_evaluated = zeros(N₋(sims.model)+N₊(sims.model))
    c = 0
    n_sim = length(sims.t)
    for i in 1:n_phases(sims.model)
        if DiscretisedFluidQueues._has_left_boundary(sims.model.S,i)
            c += 1
            pm_evaluated[c] = sum(sims.X[sims.φ.==i]==0.0)/n_sim
        end
    end
    c = 0
    for i in 1:n_phases(sims.model)
        if DiscretisedFluidQueues._has_right_boundary(model.S,i)
            c += 1
            pm_evaluated[c] = sum(sims.X[boot_sample.φ.==i]==sims.model.b)/n_sim
        end
    end
    return pm_evaluated
end
function _eval_sim(boot_sample,distribution_function,grid)
    evaluated = distribution_function(boot_sample).(grid,(1:n_phases(boot_sample.model))')
    pm_evaluated = sim_point_masses(sims)
    return evaluated, pm_evaluated
end

function compute_errors_bootstrap_ci(models,sims::Simulation,distribution_function_model,distribution_function_sim,err_fun,n_boot,grid,_phases,rng,p)
    orders = sort([l for l in keys(models)])
    k = keys(models[orders[1]])
    err_data = Matrix{Vector{Float64}}(undef,length(orders),length(k))
    n_sims = length(sims.t)
    boots = zeros(n_boot,length(orders),length(k))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        evaluated, pm_evaluated = _eval_sim(boot_sample,distribution_function_sim,grid)
        for (c_o,o) in enumerate(orders)
            for (c_m,m) in enumerate(sort([i for i in k]))
                boots[n,c_o,c_m] = log10(err_fun(distribution_function_model(models[o][m]["coeffs_mapped"]),evaluated,pm_evaluated,grid,_phases))
            end
        end
    end

    for c_o in 1:length(orders)
        for c_m in 1:length(k)
            err_data[c_o,c_m] = quantile!(boots[:,c_o,c_m],p)
        end
    end

    return err_data
end 

function analysis(err_funs,models,sims,pth,n_boot,grid,_phases,titles,rng)
    df_names = sort([k[1:end-1] for k in keys(models[minimum(keys(models))])])
    if (sum(df_names.=="dg")>1)
        df_names[findlast(df_names.=="dg")] = "dg_limiter"
    end

    err_names = string.(err_funs[i][1] for i in 1:length(err_funs))
    for (c,err_fun) in enumerate(err_funs)
        err_data = compute_errors(models,sims,err_fun[1],err_fun[2],err_fun[3],grid,_phases)
        err_data = DataFrame(err_data,df_names)
        write_errs(err_data,pth,err_names[c])
        errs_ci = compute_errors_bootstrap_ci(models,sims,err_fun[2],err_fun[3],err_fun[1],n_boot,grid,_phases,rng,[0.05;0.95])
        errs_ci = DataFrame(errs_ci,df_names)
        write_errs(errs_ci,pth,err_names[c]*"ci")
        make_err_plot(models,err_data,errs_ci,pth*"/figs/"*err_names[c],titles[c])
    end
    return nothing 
end

# function compute_errors(orders,models,truth,n_evals)
#     L1_cdf_errors = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     ks_errors = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L1_pdf_errors = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L2_pdf_errors = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))

#     k = keys(models[minimum(keys(models))])

#     for (c_o,o) in enumerate(orders)
#         L1_cdf_error_row = Vector{Union{Missing, Float64}}(undef,length(k))
#         ks_error_row = Vector{Union{Missing, Float64}}(undef,length(k))
#         L1_pdf_error_row = Vector{Union{Missing, Float64}}(undef,length(k))
#         L2_pdf_error_row = Vector{Union{Missing, Float64}}(undef,length(k))
#         for (c_m,m) in enumerate(sort([i for i in k]))
#             L1_cdf_error_row[c_m] = log10(DiscretisedFluidQueues.Lp(x->models[o][m]["cdf"](x,1),x->truth.cdf(x,1),range(-eps(),10+eps(),length=n_evals)) + 
#                                                             DiscretisedFluidQueues.Lp(x->models[o][m]["cdf"](x,2),x->truth.cdf(x,2),range(-eps(),10+eps(),length=n_evals)))
#             ks_error_row[c_m] = log10(DiscretisedFluidQueues.kolmogorov_smirnov(x->models[o][m]["cdf"](x,1),x->truth.cdf(x,1),range(-eps(),10+eps(),length=n_evals)) + 
#                                                             DiscretisedFluidQueues.kolmogorov_smirnov(x->models[o][m]["cdf"](x,2),x->truth.cdf(x,2),range(-eps(),10+eps(),length=n_evals)))
#             L1_pdf_error_row[c_m] = log10(DiscretisedFluidQueues.Lp(x->models[o][m]["pdf"](x,1),x->truth.pdf(x,1),range(-eps(),10+eps(),length=n_evals)) + 
#                                                             DiscretisedFluidQueues.Lp(x->models[o][m]["pdf"](x,2),x->truth.pdf(x,2),range(-eps(),10+eps(),length=n_evals)))
#             L2_pdf_error_row[c_m] = log10(DiscretisedFluidQueues.Lp(x->models[o][m]["pdf"](x,1),x->truth.pdf(x,1),range(-eps(),10+eps(),length=n_evals),2) + 
#                                                             DiscretisedFluidQueues.Lp(x->models[o][m]["pdf"](x,2),x->truth.pdf(x,2),range(-eps(),10+eps(),length=n_evals),2))
#         end
#         push!(L1_cdf_errors,L1_cdf_error_row)
#         push!(ks_errors,ks_error_row)
#         push!(L1_pdf_errors,L1_pdf_error_row)
#         push!(L2_pdf_errors,L2_pdf_error_row)
#     end
#     for e in (L1_cdf_errors,ks_errors,L1_pdf_errors,L2_pdf_errors)
#         for i in 1:size(e,1), j in size(e,2)
#             if ismissing(e[i,j])||isnan(e[i,j])
#                 e[i,j]=missing
#             end
#         end
#     end
#     return (L1_cdf_errors,ks_errors,L1_pdf_errors,L2_pdf_errors)
# end

function write_errs(err,pth,err_names)
    file = pth*"/data/"*string(err_names)
    CSV.write(file*".csv",err)
    return nothing
end

std_plot!(p,fun,args...; kwargs...) = plot!(
        p,
        fun,
        _xlims[1], _xlims[2];
        # ylim=(-0.1,2.1),
        xticks=false,
        yticks=false,
        grid=false,
        label=false,
        kwargs...,
    )

function make_df_plots!(models,map_to_fun,truth,truth_ci,pth,_xlims,_ylims,_xticks,_yticks,_phases)
    orders = sort([l for l in keys(models)])
    k = keys(models[orders[1]])
    p1 = plot(layout = (length(k),length(orders))) 
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(sort([i for i in k]))
            for i in _phases
                x_vals = range(_xlims[1],_xlims[2],length=200)
                y_vals = truth.(x_vals,i)
                ci = truth_ci.(x_vals,i) 
                lwr = y_vals-[ci[n][1] for n in 1:length(ci)]
                upr = [ci[n][2] for n in 1:length(ci)]-y_vals
                plot!(p1.layout.grid[c_m,c_o],x_vals,y_vals,ribbon=(lwr,upr),color=i+maximum(_phases),label=false,alpha=0.5)
                y_vals = map_to_fun(models[o][m]["coeffs_mapped"]).(x_vals,i)
                plot!(p1.layout.grid[c_m,c_o],x_vals,y_vals,color=i,label=false,alpha=0.5)
            end
            if c_o==1
                plot!(p1.layout.grid[c_m,c_o],ylabel=uppercase(string(m)[1:end-1]))
                plot!(p1.layout.grid[c_m,c_o],yticks=_yticks)
            end
            if c_m==1
                plot!(p1.layout.grid[c_m,c_o],title=string(o))
            end
            if c_m==length(k)
                plot!(p1.layout.grid[c_m,c_o],xticks=_xticks)
            end
            plot!(p1.layout.grid[c_m,c_o],ylims=_ylims)
        end
    end
    savefig(p1,pth*".svg")
    return nothing 
end

function make_err_plot(models,err_data,errs_ci,pth,title)
    orders = sort([l for l in keys(models)])
    q = plot()
    for (c_m,m) in enumerate(names(err_data))
        # idx = [!(ismissing(err[i,m])||ismissing(errs_ci[i,m][1])||ismissing(errs_ci[i,m][2])) for i in 1:length(orders)]
        @show lwr = [err_data[i,m]-errs_ci[i,m][1] for i in 1:length(orders)]
        @show upr = [errs_ci[i,m][2]-err_data[i,m] for i in 1:length(orders)]
        plot!(q,orders,err_data[:,m],ribbon=(lwr,upr),label=m,color=c_m)
    end
    plot!(q,xlabel="Order",ylabel="log₁₀ Error",title=title,legend=:outertopright)
    savefig(q,pth*".svg")
    return nothing 
end
function permute_order1_coeffs(coeffs,model,o)
    return [
        coeffs[1:N₋(model)];
        permutedims(reshape(coeffs[N₋(model)+1:end-N₊(model)],n_phases(model),o,:),[2,1,3])[:];
        coeffs[end-N₊(model)+1:end];
    ]
end
# function make_plots!(orders,models,truth,truth_ci,err_data,errs_ci,pth,_xlims,_ylims,_xticks,_yticks,_phases)
    

#     k = keys(models[minimum(keys(models))])

#     p1 = plot(layout = (length(k),length(orders))) 
#     # q1 = plot(layout = (length(k),length(orders)))

#     # what a truely terrible piece of code! 
#     for (c_o,o) in enumerate(orders)
#         for (c_m,m) in enumerate(sort([i for i in k]))
#             for i in _phases
#                 # @eval pm = sum($(Symbol(m,"_π",o))[[1:2;end-1:end]])
#                 # std_plot!(p1.layout.grid[c_m,c_o],x->models[o][m]["pdf"](x,i),color=i)
#                 x_vals = range(_xlims[1],_xlims[2],length=200)
#                 y_vals = truth.pdf.(x_vals,i)
#                 ci = truth_ci.pdf.(x_vals,i) 
#                 lwr = y_vals-[ci[n][1] for n in 1:length(ci)]
#                 upr = [ci[n][2] for n in 1:length(ci)]-y_vals
#                 plot!(p1.layout.grid[c_m,c_o],x_vals,y_vals,ribbon=(lwr,upr),color=i+2,label=false,alpha=0.5)
#                 std_plot!(p1.layout.grid[c_m,c_o],x->truth.pdf(x,i),color=i+2)
#                 # std_plot!(p1.layout.grid[c_m,c_o],x->truth_ci.pdf(x,i)[2],linestyle=:dot,color=i+2)

#                 std_plot!(q1.layout.grid[c_m,c_o],x->models[o][m]["cdf"](x,i),color=i)
#                 std_plot!(q1.layout.grid[c_m,c_o],x->truth.cdf(x,i),color=i+2)
#                 y_vals = truth.cdf.(x_vals,i)
#                 ci = truth_ci.cdf.(x_vals,i)
#                 lwr = y_vals-[ci[n][1] for n in 1:length(ci)]
#                 upr = [ci[n][2] for n in 1:length(ci)]-y_vals
#                 plot!(q1.layout.grid[c_m,c_o],x_vals,y_vals,ribbon=(lwr,upr),color=i+2,label=false,alpha=0.5)
#                 # std_plot!(q1.layout.grid[c_m,c_o],x->truth_ci.cdf(x,i)[1],linestyle=:dot,color=i+2)
#                 # std_plot!(q1.layout.grid[c_m,c_o],x->truth_ci.cdf(x,i)[2],linestyle=:dot,color=i+2)
#             end

#             if c_o==1
#                 plot!(p1.layout.grid[c_m,c_o],ylabel=uppercase(string(m)[1:end-1]))
#                 # plot!(q1.layout.grid[c_m,c_o],ylabel=uppercase(string(m)[1:end-1]))
                
#                 plot!(p1.layout.grid[c_m,c_o],yticks=_yticks[1])
#                 # plot!(q1.layout.grid[c_m,c_o],yticks=_yticks[2])
#             end
#             if c_m==1
#                 plot!(p1.layout.grid[c_m,c_o],title=string(o))
#                 # plot!(q1.layout.grid[c_m,c_o],title=string(o))
#             end
#             if c_m==4
#                 plot!(p1.layout.grid[c_m,c_o],xticks=_xticks[1])
#                 # plot!(q1.layout.grid[c_m,c_o],xticks=_xticks[2])
#             end
#             plot!(p1.layout.grid[c_m,c_o],ylims=_ylims[1])
#             # plot!(q1.layout.grid[c_m,c_o],ylims=_ylims[2])
#         end
#     end
#     file = pth*"/figs"
#     savefig(p1,file*"/pdfs.svg")
#     savefig(q1,file*"/cdfs.svg")
#     titles = [
#         "L¹ error between the true CDF and approximations";
#         "KS error between the true CDF and approximations";
#         "L¹ error between the true PDF and approximations";
#         "L² error between the true PDF and approximations";
#     ]
#     err_names = (:L1_cdf_errors,:ks_errors,:L1_pdf_errors,:L2_pdf_errors)
#     for (c_err,err) in enumerate(err_data)
#         q = plot()
#         for (c_m,m) in enumerate(names(err))
#             idx = [!(ismissing(err[i,m])||ismissing(errs_ci[1][c_err][i,m])||ismissing(errs_ci[2][c_err][i,m])) for i in 1:length(orders)]
#             lwr = err[idx,m]-errs_ci[1][c_err][idx,m]
#             upr = errs_ci[2][c_err][idx,m]-err[idx,m]
#             plot!(q,orders[idx],err[idx,m],ribbon=(lwr,upr),label=m,color=c_m)
#             # idx = [!ismissing(errs_ci[1][c_err][i,m]) for i in 1:length(orders)]
#             # plot!(q,orders[idx],errs_ci[1][c_err][idx,m],label=false,linestyle=:dot,color=c_m)
#             # idx = [!ismissing(errs_ci[2][c_err][i,m]) for i in 1:length(orders)]
#             # plot!(q,orders[idx],errs_ci[2][c_err][idx,m],label=false,linestyle=:dot,color=c_m)
#         end
#         plot!(q,xlabel="Order", ylabel="log₁₀ Error", 
#             title=titles[c_err],
#             legend=:outertopright)
#         #display(q)
#         # #display((@__DIR__)*)
#         file = pth*"/figs/"*string(err_names[c_err])
#         savefig(q,file*".svg")
#     end 
#     return p1, q1 
# end

function cdf_bootstrap_ci(sims,n_boot,grid,_phases,p,rng)
    n_sims = length(sims.φ)
    cdfs = Array{Float64,3}(undef,n_boot,length(grid),length(_phases))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        cdfs[n,:,:] = cdf(boot_sample).(grid,_phases)
    end
    ci_lwr = Array{Float64,2}(undef,size(grid)...)
    ci_upr = Array{Float64,2}(undef,size(grid)...)
    for g in 1:length(grid)
        for i in 1:length(phases)
            ci_lwr[g,i] = quantile!(cdfs[:,g,i])
            ci_upr[g,i] = quantile!(cdfs[:,g,i])
        end
    end
    return ci_lwr, ci_upr
end

# function bootstrap_errs_ci(sims,orders,models,n_boot,n_evals,p,rng)
#     n_sims = length(sims.φ)
#     boot_errors = Array{Any,1}(undef,n_boot)
#     for n in 1:n_boot
#         idx = rand(rng,1:n_sims,n_sims)
#         boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
#         truth =  (pdf=(x,i)->missing, cdf=cdf(boot_sample))
#         boot_errors[n] = compute_errors(orders,models,truth,n_evals)
#     end
#     k = keys(models[minimum(keys(models))])

#     L1_cdf_errors_lwr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     ks_errors_lwr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L1_pdf_errors_lwr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L2_pdf_errors_lwr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L1_cdf_errors_upr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     ks_errors_upr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L1_pdf_errors_upr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     L2_pdf_errors_upr = DataFrame(DG = Vector{Union{Missing, Float64}}(undef,0), DG_Limit = Vector{Union{Missing, Float64}}(undef,0), Order_1 = Vector{Union{Missing, Float64}}(undef,0), QBDRAP = Vector{Union{Missing, Float64}}(undef,0))
#     for (c_o,o) in enumerate(orders)
#         L1_cdf_error_row_lwr = Vector{Union{Missing, Float64}}(undef,length(k))
#         ks_error_row_lwr = Vector{Union{Missing, Float64}}(undef,length(k))
#         L1_pdf_error_row_lwr = Vector{Union{Missing, Float64}}(undef,length(k))
#         L2_pdf_error_row_lwr = Vector{Union{Missing, Float64}}(undef,length(k))
#         L1_cdf_error_row_upr = Vector{Union{Missing, Float64}}(undef,length(k))
#         ks_error_row_upr = Vector{Union{Missing, Float64}}(undef,length(k))
#         L1_pdf_error_row_upr = Vector{Union{Missing, Float64}}(undef,length(k))
#         L2_pdf_error_row_upr = Vector{Union{Missing, Float64}}(undef,length(k))
#         for (c_m,m) in enumerate(sort([i for i in k]))
#             L1_cdf_error_row_lwr[c_m] = quantile(skipmissing([boot_errors[n][1][c_o,c_m] for n in 1:n_boot]),p[1])
#             ks_error_row_lwr[c_m] = quantile(skipmissing([boot_errors[n][2][c_o,c_m] for n in 1:n_boot]),p[1])
#             L1_pdf_error_row_lwr[c_m] = missing#quantile(skipmissing([boot_errors[n][3][c_o,c_m] for n in 1:n_boot]),p[1])
#             L2_pdf_error_row_lwr[c_m] = missing#quantile(skipmissing([boot_errors[n][4][c_o,c_m] for n in 1:n_boot]),p[1])
#             L1_cdf_error_row_upr[c_m] = quantile(skipmissing([boot_errors[n][1][c_o,c_m] for n in 1:n_boot]),p[2])
#             ks_error_row_upr[c_m] = quantile(skipmissing([boot_errors[n][2][c_o,c_m] for n in 1:n_boot]),p[2])
#             L1_pdf_error_row_upr[c_m] = missing#quantile(skipmissing([boot_errors[n][3][c_o,c_m] for n in 1:n_boot]),p[2])
#             L2_pdf_error_row_upr[c_m] = missing#quantile(skipmissing([boot_errors[n][4][c_o,c_m] for n in 1:n_boot]),p[2])
#         end
#         push!(L1_cdf_errors_lwr,L1_cdf_error_row_lwr)
#         push!(ks_errors_lwr,ks_error_row_lwr)
#         push!(L1_pdf_errors_lwr,L1_pdf_error_row_lwr)
#         push!(L2_pdf_errors_lwr,L2_pdf_error_row_lwr)
#         push!(L1_cdf_errors_upr,L1_cdf_error_row_upr)
#         push!(ks_errors_upr,ks_error_row_upr)
#         push!(L1_pdf_errors_upr,L1_pdf_error_row_upr)
#         push!(L2_pdf_errors_upr,L2_pdf_error_row_upr)
#     end

#     return (L1_cdf_errors_lwr,ks_errors_lwr,L1_pdf_errors_lwr,L2_pdf_errors_lwr), (L1_cdf_errors_upr,ks_errors_upr,L1_pdf_errors_upr,L2_pdf_errors_upr)
# end