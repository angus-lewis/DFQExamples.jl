include("default_params.jl")
include("reflecting_model/model_def.jl")

include("../helper_functions.jl")

# initial condition
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

function map_to_t(d0,B,approx_type,t,h=0.005,limit=false,kwargs...)
    if (approx_type==:limiter)
        dt = integrate_time(d0,B,t,StableRK4(h); limiter=GeneralisedMUSCL)
    end
    return dt
end

function make_approx_data()
    pth = mkpath((@__DIR__)*"/reflecting_model/data/exp/data")
    mkpath(pth*"/coeffs")
    mkpath(pth*"/cdfs")

    h0 = 1
    for k in 1:2:21
        h = h0/((k+1)/2)
        o = 2
        nodes = 0:h:10

        mesh = DGMesh(nodes,o)
        
        dq = DiscretisedFluidQueue(model,mesh)

        d0 = d0_map_exp(dq)

        B = build_full_generator(dq)

        dt = map_to_t(d0,B.B,:limiter,transient_time)

        CSV.write(
            pth*"/coeffs/order_"*string(k)*"_model_limiter.csv",
            DataFrame(a=dt.coeffs),
        )
        x = range(0,model.b,length=n_err_evals)
        phase_1_evaluated = cdf(dt).(x,1)
        phase_2_evaluated = cdf(dt).(x,2)
        CSV.write(
            pth*"/cdfs/order_"*string(k)*"_model_limiter.csv",
            DataFrame(x=x,
                phase_1 = phase_1_evaluated,
                phase_2 = phase_2_evaluated,
            ),
        )
    end
    return nothing
end

make_approx_data()