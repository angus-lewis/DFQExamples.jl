include((@__DIR__)*"/../preamble.jl")

function map_to_t(d0,B,approx_type,t,h=0.005,limit=false,kwargs...)
    if (approx_type==:dg)&&(limit)
        dt = integrate_time(d0,B,t,StableRK4(h); limiter=GeneralisedMUSCL)
    elseif (approx_type==:dg)&&(!limit)
        dt = integrate_time(d0,B,t,StableRK4(h); limiter=NoLimiter)
    else 
        dt = integrate_time(d0,B,t,StableRK4(h))
    end
    return dt
end

function run_transient(m_str,ic_str,ic,d0_map,rng)
    # load model and parameters 
    model_str = (@__DIR__)*"/"*m_str
    include(model_str*"/model_def.jl")
    include((@__DIR__)*"/default_params.jl")

    # make folders to write to 
    pth = model_str*"/transient_distribution/"*ic_str
    mkpath(pth*"/data")
    mkpath(pth*"/figs")

    # construct all the approximations
    models = make_approximations(orders,approx_types,d0_map,map_to_t)
    sims = simulate_model(model,fixed_time(transient_time),n_sims,ic,pth,rng)

    # evaluate models 
    _phases = 1:2
    grid = range(0.0,sims.model.b,length=n_err_evals)
    nodes = models[orders[1]]["dg1"]["dq"].mesh.nodes
    centres = (nodes[1:end-1] + nodes[2:end])./2
    @evaluate!(models,(cdf,pdf,cell_probs),gird,_phases)
    # @eval models=$models
    # evlauate sims
    sims_evaluated = @evaluate(sims,(cdf,),grid,_phases)
    sims_evaluated["cell_probs"] = cell_probs(sims,nodes).(centres,_phases)
    # @eval sims=$sims
    truth_lwr, truth_upr = cdf_bootstrap_ci(sims,n_boot,grid,_phases,(0.05,0.95),rng)
    truth = cdf(sims)

    # parameters for plotting
    _ylims = (-0.01,0.65)#(-0.01,0.4)
    _yticks = 0:0.1:0.65#(0:0.1:0.4,
    _xlims = (-eps(),model.b+eps())
    _xticks = 0:5:10
    make_df_plots!(models,cdf,truth,cdf_ci,pth*"/figs/cdfs",_xlims,_ylims,_xticks,_yticks,_phases)
    _ylims = (-0.01,0.4)
    _yticks = 0:0.1:0.4
    make_df_plots!(models,pdf,(x,i)->NaN,(x,i)->(NaN,NaN),pth*"/figs/pdfs",_xlims,_ylims,_xticks,_yticks,_phases)

    # analysis on cdfs and cell_probs only 
    # @eval DiscretisedFluidQueues.cdf(s::Simulation,v) = DiscretisedFluidQueues.cdf(s)
    err_funs = ((kolmogorov_smirnov,identity,cdf),(L1,cdf,cdf),(L2,cdf,cdf)) 
    titles = [
        "KS error between the true CDF and approximations";
        "L¹ error between the true CDF and approximations";
        "L² error between the true CDF and approximations";
    ]
    analysis(err_funs,models,sims,pth,n_boot,grid,_phases,titles,rng)

    err_funs = ((L1_cell_probs,identity,x->cell_probs(x,nodes)),) 
    titles = ["L¹ error between the cell masses"]
    for o in keys(models)
        models[o]["order13"]["coeffs_mapped"] = SFMDistribution(
            permute_order1_coeffs(models[o]["order13"]["coeffs_mapped"].coeffs,model,o),
            models[o]["dg1"]["coeffs_mapped"].dq)
    end
    # @eval models=$models 
    analysis(err_funs,models,sims,pth,n_boot,grid,_phases,titles,rng)
    return nothing 
end

point_mass_initial_condition(n,rng) = (X = fill(eps(),n), φ = fill(1,n))
d0_map_point_mass(dq) = interior_point_mass(eps(),1,dq)


function exp_initial_condition(n,rng) 
    i = rand(rng,1:2,n) 
    X = Vector{Float64}(undef,n)
    for k in 1:n
        x = Inf
        while x>=model.b
            x = -log(rand(rng))
        end
        X[k]=x
    end
    return (X=X,φ=i)
end
d0_map_exp(dq) = SFMDistribution((x,i)->(i∈(1:2))*exp(-x)/(1-exp(-model.b))/2,dq)