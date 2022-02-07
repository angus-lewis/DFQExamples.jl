include((@__DIR__)*"/../preamble.jl") 
include("reflecting_model/model_def.jl")

include("default_params.jl")

# look at time to exit of the interval [0,1]
# use the bounded fluid queue, but change the boundary to be at 1

dg_model(order) = build_discretised_model(DGMesh{t},model,1.0,order)
order1_model(order) = build_discretised_model(DGMesh{t},model,1.0/order,1)
qbdrap_model(order) = build_discretised_model(FRAPMesh{t},model,1.0,order)

d0_map_point_mass(dq) = (interior_point_mass(eps(),1,dq).coeffs)

# change d0_map here to do the make_approximations for ohter initial confitions
# remember to change write directory below too
models = make_approximations(orders,approx_types,d0_map_point_mass,(args...)->args[1])

h = 0.005
pth = mkpath((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2")
pth_coeff = mkpath((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2/approx_coeffs")
pth_cdf = mkpath((@__DIR__)*"/reflecting_model/transient_distribution/point_mass/datat2/approx_cdf")
x_vals = range(0,10;length=n_err_evals)
for k1 in keys(models)
    for k2 in keys(models[k1])
        models[k1][k2]["coeffs_mapped"] = Vector{Float64}(JSON.parsefile(
            (@__DIR__)*
            "/reflecting_model/transient_distribution/point_mass/data/models/order_"*string(k1)*
            "/model_"*k2*"/coeffs_mapped.json",
        ))
        @show k1
        @show k2
        if k2=="dg2"
            models[k1][k2]["coeffs_mapped"] = integrate_time(SFMDistribution(models[k1][k2]["coeffs_mapped"],models[k1][k2]["dq"]),models[k1][k2]["B"],0.1,StableRK4(h); limiter=GeneralisedMUSCL)
        elseif k2=="dg1"
            models[k1][k2]["coeffs_mapped"] = integrate_time(SFMDistribution(models[k1][k2]["coeffs_mapped"],models[k1][k2]["dq"]),models[k1][k2]["B"],0.1,StableRK4(h); limiter=NoLimiter)
        else 
            models[k1][k2]["coeffs_mapped"] = integrate_time(SFMDistribution(models[k1][k2]["coeffs_mapped"],models[k1][k2]["dq"]),models[k1][k2]["B"],0.1,StableRK4(h))
        end
        CSV.write(pth_coeff*"/order_"*string(k1)*"model_"*k2*".csv",DataFrame(coeffs=models[k1][k2]["coeffs_mapped"]))
        cdf_evaluated = DataFrame(
            x=x_vals,
            phase_1 = cdf(models[k1][k2]["coeffs_mapped"]).(x,1),
            phase_2 = cdf(models[k1][k2]["coeffs_mapped"]).(x,2),
        )
        CSV.write(pth_cdf*"/order_"*string(k1)*"model_"*k2*".csv",cdf_evaluated)
    end
end




