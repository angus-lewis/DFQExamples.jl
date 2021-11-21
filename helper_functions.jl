"""
    Given a Mesh, a BoundedFluidQueue, a grid spacing hx, and order construct 
a DiscretisedFluidQueue as specified by Mesh, and a Generator. By default 
also constructs an initial condition which approximates a point mass at x=eps(),
in phase i=1
"""
function build_discretised_model(type::Type{<:Mesh{T}},model,hx,order) where T 
    nodes = 0.0:hx:model.b
    mesh = type(nodes,order)
    dq = DiscretisedFluidQueue(model,mesh)
    B = build_full_generator(dq)
    d0 = interior_point_mass(eps(),1,dq)
    return d0, B, dq
end

t = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
# some functions to help buid approximations en masse
dg_model(order) = build_discretised_model(DGMesh{t},model,1.0,order)
order1_model(order) = build_discretised_model(DGMesh{t},model,1.0/order,1)
qbdrap_model(order) = build_discretised_model(FRAPMesh{t},model,1.0,order)

"""
Build and save all of the approximations as specified by `orders` and `approx_types`.
`orders` is a vector of odd integers. `approx_types` is a vector of symbols. If `:dg`
appears twice in `approx_types` then a flux limiter is applied to the second one. 

`d0_map` is a function which take a `DiscretisedFluidQueue` and maps it to an initial 
condition. 

`a_map` is a mapping to apply to the approximation which returns the information of interest. 
For example, `map_to_t(d0,B,approx_type,t,h=0.005,limit=false,kwargs...)`.

Construct a Dict of `models` with keys `orders`. `models[o]` is another Dict
with keys `approx_types` and values as follows: 
    "pdf"=>the_pdf,
    "cdf"=>the_cdf,
    "coeffs_mapped"=>dt,
    "coeffs_0"=>d0,
    "dq"=>dq,
    "B"=>B,
"""
function make_approximations(orders,approx_types,d0_map,a_map=identity)
    models = Dict{Any,Dict{Any,Dict}}()
    for (c_o,o) in enumerate(orders)
        tmp_dict = Dict()
        dg_counter = 0
        for (c_m,m) in enumerate(approx_types)

            ~,B,dq = @eval $(Symbol(m,"_model"))($o)

            d0 = d0_map(dq)

            (string(m)=="dg")&&(dg_counter+=1)
            # if two dg methods are specified, use a limiter with the second one
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

"""
Automatically save the output of `evalate!` under the function name (as a string) as the key.
"""
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

"""
Evaluate the functions in `funs` for all of the `models` at the x-values in `grid`
and phases in `_phases`. 

If `f` is in `funs` then we evaluate `f(models[o][m]).(grid[:],_phases[:]')` where 
`m` and `o` are keys specifying the approximation type and order, respectively. 
The output is stored in `models[o][m][fun_name]`.
"""
function evaluate!(models,funs::Tuple,grid::AbstractVector,
    _phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple})

    orders = sort([l for l in keys(models)])
    k = sort([k for k in keys(first(models)[2])])

    # all the models are the same, so just grab one of them
    model = models[orders[1]][k[1]]["coeffs_mapped"].dq.model

    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(k)
            for (c_f,f) in enumerate(funs)
                temp_fun = f(models[o][m]["coeffs_mapped"])
                models[o][m][fun_names[c_f]*"_evaluated"] = temp_fun.(grid[:],_phases[:]')
                models[o][m]["boundary_masses"] = models[o][m]["coeffs_mapped"][[1:N₋(model);end-N₊(model)+1:end]]
            end
        end
    end
    return nothing
end


"""
Evaluate `evaluate(models,funs,grid,phases)` where `phases` are all the phases in the model
"""
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
"""
Evaluate `evaluate(models,funs,grid,phases)` where `phases` are all the phases in the model
"""
function evaluate!(models,funs::Tuple,grid::AbstractVector,fun_names::Union{AbstractVector{String},Tuple}) 
    o = sort([l for l in keys(models)])
    k = keys(models[o[1]])
    model = models[o[1]][k[1]]["coeffs_mapped"].dq.model
    evaluate!(models,funs,grid,1:n_phases(model),fun_names)
    return nothing
end
"""
Call evaluate(sims::Simulation,funs::Tuple,grid::AbstractVector,_phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple})
where `fun_names` are the names of the functions supplied in `funs`
"""
macro evaluate(sims,funs,grid,_phases)
    # esc(funs).args[1].args -- esc() returns an expression, the first .args gets the args of the output of esc()
    # the second .args return the function names in `funs` as a vector of symbols
    return quote 
        evaluate($(esc(sims)),$(esc(funs)),$(esc(grid)),$(esc(_phases)),$(string.(esc(funs).args[1].args)))
    end
end
"""
Evaluate the functions in `funs` for the `sim` at the x-values in `grid`
and phases in `_phases`. 

If `f` is in `funs` then we evaluate `f(models[o][m]).(grid[:],_phases[:]')` where 
`m` and `o` are keys specifying the approximation type and order, respectively. 
The output is stored in `models[o][m][fun_name]`.
"""
function evaluate(sims::Simulation,funs::Tuple,grid::AbstractVector,_phases::AbstractVector{Int},fun_names::Union{AbstractVector{String},Tuple})
    evaluated = Dict{String,Any}()
    for (c_f,f) in enumerate(funs)
        temp_fun = f(sims)
        evaluated[fun_names[c_f]*"_evaluated"] = temp_fun.(grid[:],_phases[:]')
        evaluated["boundary_masses"] = sim_point_masses(sims)
    end
    return evaluated
end
macro evaluate(sims,funs,grid)
    return quote 
        evaluate($(esc(sims)),$(esc(funs)),$(esc(grid)),$(string.(esc(funs).args[1].args)))
    end
end
evaluate(sims::Simulation,funs::Tuple,grid::AbstractVector,fun_names::Union{AbstractVector{String},Tuple}) =
    evaluate(sims,funs,grid,1:n_phases(sims.model),fun_names)

"""
evaluate `cell_probs(model).(centres[:],_phases[:]')` for all models as specified in the Dict 
`models` (as outpud from `make_approximations()`)
"""
function evaluate_cell_probs!(models,nodes,_phases)
    dict_names = sort([k for k in keys(first(models)[2])])
    orders = sort([l for l in keys(models)])
    centres = (nodes[1:end-1]+nodes[2:end])./2
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(dict_names)
            if m[1:end-1]!="order1"
                models[o][m]["cell_probs_evaluated"] = cell_probs(models[o][m]["coeffs_mapped"]).(centres[:],_phases[:]')
            else
                # need to map the coefficients of the order1 method to have the same cell_width
                # as the dg method.
                d = SFMDistribution(
                    permute_order1_coeffs(models[o][m]["coeffs_mapped"].coeffs,
                        models[o][m]["coeffs_mapped"].dq.model,
                        o,
                    ),
                    models[o]["dg1"]["coeffs_mapped"].dq
                )
                models[o][m]["cell_probs_evaluated"] = cell_probs(d).(centres[:],_phases[:]')
            end
        end
    end
    return nothing 
end 

"""
Simulate, `n_sim` times the `BoundedFluidQueue` as specified by `model` until the `stopping_time` with 
initial points specified by `init` where rows of init are samples of [X0, φ0]
"""
function simulate_model(model,stopping_time,n_sims,init,pth,rng::Random.AbstractRNG=Randon.default_rng())
    SFM0 = init(n_sims,rng)
    sims = simulate(model,stopping_time,SFM0,rng)
    df = DataFrame(t=sims.t,φ=sims.φ,X=sims.X,n=sims.n)
    CSV.write(pth*"/data/sim.csv",df)
    return sims
end

"""
compute the values of the point masses as approximated by the data in sims.
"""
function sim_point_masses(sims::Simulation)
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
            pm_evaluated[c] = sum(sims.X[sims.φ.==i]==sims.model.b)/n_sim
        end
    end
    return pm_evaluated
end

function bootstrap_error_ci(f1_evaluated,f1_pm,boot_samples_evaluated,boot_samples_pm_evaluated,grid,error_metric,p)
    boot_errors = Vector{Float64}(undef,size(boot_samples_evaluated,1))
    for n in 1:size(boot_samples_evaluated,1)
        boot_errors[n] = error_metric(f1_evaluated,f1_pm,boot_samples_evaluated[n,:,:],boot_samples_pm_evaluated[n,:,:],grid)
    end
    return (quantile!(boot_errors,p)..., boot_errors)
end

function analysis(err_funs,models,truth,grid,p,titles,pth; w=true)
    dict_names = sort([k for k in keys(first(models)[2])])
    df_names = copy(dict_names)
    if (sum(df_names.=="dg")>1)
        df_names[findlast(df_names.=="dg")] = "dg_limiter"
    end
    err_names = string.(err_funs[i][1] for i in 1:length(err_funs))
    orders = sort([k for k in keys(models)])
    for (c_f,f) in enumerate(err_funs)
        error_metric = f[1]
        apply_to = f[2]
        err_data = Matrix{Float64}(undef,length(orders),length(dict_names))
        errs_ci_lwr = Matrix{Float64}(undef,length(orders),length(dict_names))
        errs_ci_upr = Matrix{Float64}(undef,length(orders),length(dict_names))
        (w)&&(specific_path=mkpath(pth*"/data/bootstrap_errors/"*err_names[c_f]*"/"*apply_to*"/"))
        for (c_o,o) in enumerate(orders)
            for (c_m,m) in enumerate(dict_names)
                err_data[c_o,c_m] = error_metric(models[o][m][apply_to*"_evaluated"],models[o][m]["boundary_masses"],truth[apply_to*"_evaluated"],truth["boundary_masses"],grid)
                errs_ci_lwr[c_o,c_m], errs_ci_upr[c_o,c_m], bootstrap_data = 
                    bootstrap_error_ci(
                        models[o][m][apply_to*"_evaluated"],models[o][m]["boundary_masses"],
                        truth["bootstrap_samples"][apply_to],truth["bootstrap_samples"]["boundary_masses"],
                        grid,error_metric,p,
                    )#error_metric(models[o][m][apply_to*"_evaluated"],models[o][m]["boundary_masses"],truth[apply_to*"_ci_lwr"],truth["boundary_masses"*"_ci_lwr"],grid)
                # errs_ci_upr[c_o,c_m]  = error_metric(models[o][m][apply_to*"_evaluated"],models[o][m]["boundary_masses"],truth[apply_to*"_ci_upr"],truth["boundary_masses"*"_ci_upr"],grid)
                (w)&&JSON.write(specific_path*"/order_"*string(o)*"_model_"*m*".json",JSON.json(bootstrap_data))
            end
        end
        err_data = DataFrame(err_data,df_names)
        write_errs(err_data,pth,err_names[c_f])
        errs_ci_lwr = DataFrame(errs_ci_lwr,df_names)
        write_errs(errs_ci_lwr,pth,err_names[c_f]*"_ci_lwr")
        errs_ci_upr = DataFrame(errs_ci_upr,df_names)
        write_errs(errs_ci_upr,pth,err_names[c_f]*"ci_upr")
        (w)&&make_err_plot(models,err_data,errs_ci_lwr,errs_ci_upr,pth*"/figs/"*err_names[c_f],titles[c_f])
    end
    return nothing 
end

function make_err_plot(models,err_data,errs_ci_lwr,errs_ci_upr,pth,title)
    orders = sort([l for l in keys(models)])
    q = plot()
    for (c_m,m) in enumerate(names(err_data))
        lwr = err_data[:,m]-errs_ci_lwr[:,m]
        upr = errs_ci_upr[:,m]-err_data[:,m]
        plot!(q,orders,err_data[:,m],ribbon=(lwr,upr),label=m,color=c_m)
    end
    plot!(q,xlabel="Order",ylabel="log₁₀ Error",title=title,legend=:outertopright)
    savefig(q,pth*".svg")
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

function make_df_plots!(models,truth,map_name,grid,_phases,pth,_xlims,_ylims,_xticks,_yticks)
    orders = sort([l for l in keys(models)])
    k = keys(models[orders[1]])
    p1 = plot(layout = (length(k),length(orders))) 
    t = keys(truth)
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(sort([i for i in k]))
            for (c_i,i) in enumerate(_phases)
                (map_name*"_evaluated"∈t)&&(y_vals = truth[map_name*"_evaluated"][:,c_i])
                if (map_name*"_ci_lwr")∈t
                    lwr = y_vals-truth[map_name*"_ci_lwr"][:,c_i]
                    upr = truth[map_name*"_ci_upr"][:,c_i]-y_vals
                    plot!(p1.layout.grid[c_m,c_o],grid,y_vals,ribbon=(lwr,upr),color=i+maximum(_phases),label=false,xticks=false,alpha=0.5)
                elseif (map_name∈t)
                    plot!(p1.layout.grid[c_m,c_o],grid,y_vals,color=i+maximum(_phases),label=false,xticks=false,alpha=0.5)
                end
                y_vals = models[o][m][map_name*"_evaluated"][:,c_i]
                plot!(p1.layout.grid[c_m,c_o],grid,y_vals,color=i,xticks=false,label=false,alpha=0.5)
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
function ci_from_samples(samples,p)
    ci_lwr = Array{Float64,2}(undef,size(samples)[2:end]...)
    ci_upr = Array{Float64,2}(undef,size(samples)[2:end]...)
    for g in 1:size(samples,2)
        for i in 1:size(samples,3)
            ci_lwr[g,i], ci_upr[g,i] = quantile!(samples[:,g,i],p)
        end
    end
    return ci_lwr, ci_upr
end
function cdf_bootstrap_ci(sims,n_boot,grid,_phases,p,rng)
    n_sims = length(sims.φ)
    cdfs = Array{Float64,3}(undef,n_boot,length(grid),length(_phases))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        cdfs[n,:,:] = cdf(boot_sample).(grid[:],_phases[:]')
    end
    ci_lwr, ci_upr = ci_from_samples(cdfs,p)
    return ci_lwr, ci_upr
end
function cell_probs_bootstrap_ci(sims,n_boot,nodes,_phases,p,rng)
    n_sims = length(sims.φ)
    centres = (nodes[1:end-1]+nodes[2:end])/2.0
    cell_probs_evald = Array{Float64,3}(undef,n_boot,length(centres),length(_phases))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        cell_probs_evald[n,:,:] = cell_probs(boot_sample,nodes).(centres[:],_phases[:]')
    end
    ci_lwr, ci_upr = ci_from_samples(cell_probs_evald,p)
    return ci_lwr, ci_upr
end
function point_mass_boostrap_ci(sims,n_boot,p,rng)
    n_sims = length(sims.φ)
    pm_samps = Array{Float64,3}(undef,n_boot,1,N₋(sims.model)+N₊(sims.model))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        pm_samps[n,1,:] = sim_point_masses(boot_sample)
    end
    ci_lwr, ci_upr = ci_from_samples(pm_samps,p)
    return ci_lwr, ci_upr
end

function bootstrap_cis(sims,n_boot,grid,nodes,_phases,p,rng)
    n_sims = length(sims.φ)
    cdfs = Array{Float64,3}(undef,n_boot,length(grid),length(_phases))
    centres = (nodes[1:end-1]+nodes[2:end])/2.0
    cell_probs_evald = Array{Float64,3}(undef,n_boot,length(centres),length(_phases))
    pm_samps = Array{Float64,3}(undef,n_boot,1,N₋(sims.model)+N₊(sims.model))
    for n in 1:n_boot
        idx = rand(rng,1:n_sims,n_sims)
        boot_sample = Simulation(sims.t[idx],sims.φ[idx],sims.X[idx],sims.n[idx],sims.model)
        cdfs[n,:,:] = cdf(boot_sample).(grid[:],_phases[:]')
        cell_probs_evald[n,:,:] = cell_probs(boot_sample,nodes).(centres[:],_phases[:]')
        pm_samps[n,1,:] = sim_point_masses(boot_sample)
    end
    cdf_ci_lwr, cdf_ci_upr = ci_from_samples(cdfs,p)
    cell_probs_ci_lwr, cell_probs_ci_upr = ci_from_samples(cell_probs_evald,p)
    pm_ci_lwr, pm_ci_upr = ci_from_samples(pm_samps,p)
    return Dict{String,Any}(
        "bootstrap_samples" => Dict("cdf"=>cdfs,"cell_probs"=>cell_probs_evald,"boundary_masses"=>pm_samps),
        "cdf_ci_lwr" => cdf_ci_lwr, 
        "cdf_ci_upr" => cdf_ci_upr,
        "cell_probs_ci_lwr" => cell_probs_ci_lwr, 
        "cell_probs_ci_upr" => cell_probs_ci_upr,
        "boundary_masses_ci_lwr" => pm_ci_lwr[:], 
        "boundary_masses_ci_upr" => pm_ci_upr[:],
    )
end

function write_models(models,pth)
    orders = sort([l for l in keys(models)])
    k = sort([k for k in keys(first(models)[2])])
    for (c_o,o) in enumerate(orders)
        for (c_m,m) in enumerate(k)
            temp_path = mkpath(pth*"/data/models/order_"*string(o)*"/model_"*string(m)*"/")
            for var in keys(models[o][m])
                try 
                    JSON.write(temp_path*"/"*var,JSON.json(models[o][m][var]))
                catch 
                    JSON.write(temp_path*"/"*var,JSON.json("failed to write "*var))
                end
            end
        end
    end
    return nothing 
end
function write_truth(truth,pth)
    temp_path = mkpath(pth*"/data/sims/")
    for var in keys(truth)
        try 
            JSON.write(temp_path*"/"*var,JSON.json(truth[var]))
        catch 
            JSON.write(temp_path*"/"*var,JSON.json("failed to write "*var))
        end
    end
    return nothing 
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