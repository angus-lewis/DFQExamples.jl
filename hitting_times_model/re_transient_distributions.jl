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

function re_run_transient(
    orders,approx_types,transient_time,n_sims,n_err_evals,n_boot,qtiles,
    model,m_str,ic,ic_str,d0_map,rng,
)

    # load model and parameters 
    model_str = (@__DIR__)*"/"*m_str

    # make folders to write to 
    pth = model_str*"/transient_distribution/"*ic_str
    mkpath(pth*"/data")
    mkpath(pth*"/figs")

    # construct all the approximations
    models = make_approximations(orders,approx_types,d0_map,map_to_t)
    sims = simulate_model(model,fixed_time(transient_time),n_sims,ic,pth,rng)

    # evaluate models 
    _phases = 1:4
    grid = range(0.0,sims.model.b,length=n_err_evals)
    nodes = 0:1:10#models[orders[1]]["dg1"]["dq"].mesh.nodes
    centres = (nodes[1:end-1] + nodes[2:end])./2
    # @evaluate!(models,(cdf,pdf),grid,_phases)
    evaluate_cell_probs!(models,nodes,_phases)

    # evlauate sims
    truth = @evaluate(sims,(cdf,),grid,_phases)
    truth["cell_probs_evaluated"] = cell_probs(sims,nodes).(centres[:],_phases[:]')
    truth = merge!(truth,bootstrap_cis(sims,n_boot,grid,nodes,_phases,qtiles,rng))

    # parameters for plotting
    _ylims = (-0.01,0.65)#(-0.01,0.4)
    _yticks = 0:0.1:0.65#(0:0.1:0.4,
    _xlims = (-eps(),model.b+eps())
    _xticks = 0:5:10
    make_df_plots!(models,truth,"cdf",grid,_phases,pth*"/figs/cdfs",_xlims,_ylims,_xticks,_yticks)
    _ylims = (-0.01,0.4)
    _yticks = 0:0.1:0.4
    make_df_plots!(models,truth,"pdf",grid,_phases,pth*"/figs/pdfs",_xlims,_ylims,_xticks,_yticks)

    # analysis on cdfs and cell_probs only 
    err_funs = ((kolmogorov_smirnov,"cdf"),(L1,"cdf"),(L2,"cdf"),(L1_cell_probs,"cell_probs")) 
    titles = [
        "KS error between the true CDF and approximations";
        "L¹ error between the true CDF and approximations";
        "L² error between the true CDF and approximations";
        "L¹ error between the cell masses";
    ]
    analysis(err_funs,models,truth,grid,qtiles,titles,pth;w=true)
    # write all the data for back up
    write_models(models,pth)
    write_truth(truth,pth)
    return nothing 
end
