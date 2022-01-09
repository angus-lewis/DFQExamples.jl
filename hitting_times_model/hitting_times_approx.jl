include((@__DIR__)*"/../preamble.jl") 

include("absorbing_model/model_def.jl")
include("default_params.jl")

# look at time to exit of the interval [0,1]
# use the bounded fluid queue, but change the boundary to be at 1
model = BoundedFluidQueue(
    model.T[1:2,1:2],rates(model)[1:2],[0.0 1.0],[1.0 0.0],1.0,
)

dg_model(order) = build_discretised_model(DGMesh{t},model,1.0/3,order)
order1_model(order) = build_discretised_model(DGMesh{t},model,1.0/3/order,1)
qbdrap_model(order) = build_discretised_model(FRAPMesh{t},model,1.0/3,order)

d0_map_point_mass(dq) = (interior_point_mass(eps(),1,dq).coeffs)
d0_map_exp(dq) = (SFMDistribution((x,i)->(i∈(1:2))*exp(-x)/(1-exp(-model.b))/2,dq).coeffs)

models = make_approximations(orders,approx_types,d0_map_exp,(args...)->args[1])

for k1 in keys(models)
    for k2 in keys(models[k1])
        models[k1][k2]["B"][1,:] = spzeros(1,size(models[k1][k2]["B"],2))
        models[k1][k2]["B"][end,:] = spzeros(1,size(models[k1][k2]["B"],2))
        droptol!(models[k1][k2]["B"],eps()*1000)
    end
end

pth = mkpath((@__DIR__)*"/hitting_times/data/exp")

h = 0.005
t_vec = 0.0:h:10.0
hitting_times_cdf = zeros(length(t_vec),3) 
hitting_times_cdf[:,1] = t_vec
for k1 in keys(models)
    for k2 in keys(models[k1])
        for (c,t) in enumerate(t_vec)
            models[k1][k2]["coeffs_mapped"] = 
                DiscretisedFluidQueues._integrate(
                    models[k1][k2]["coeffs_mapped"],
                    models[k1][k2]["B"],
                    h,
                    StableRK4(h),
                    (k2!="dg2" ? identity : (x->GeneralisedMUSCL.fun(x,GeneralisedMUSCL.generate_params(models[k1][k2]["dq"])...))),
                )
            hitting_times_cdf[c,2:3] = models[k1][k2]["coeffs_mapped"][[1;end]]
        end
        CSV.write(
            pth*"/order_"*string(k1)*"_model_"*string(k2)*".csv",
            DataFrame(
                t=hitting_times_cdf[:,1],
                phase_2=hitting_times_cdf[:,2],
                phase_1=hitting_times_cdf[:,3],
            )
        )
    end
end

