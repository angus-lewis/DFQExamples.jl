include("default_params.jl")
include("model_def.jl")

include("../helper_functions.jl")

orders = orders[2:3:end]
# initial condition
d0_map_point_mass(dq) = (interior_point_mass(5.0,3,dq).coeffs)

# build the generators
models = make_approximations(orders,approx_types,d0_map_point_mass,(args...)->args[1])

# points at which to evaluate the first return cdf
x_vals = range(0.0,1.6;length=n_err_evals)

for k1 in orders
    for k2 in keys(models[k1])
        if k2!="dg2"
            @show k1
            @show k2
            println()
            # evaluate psi for all approximations
            _ffq = ffq(models[k1][k2]["dq"].mesh,Δ(models[k1][k2]["dq"],1))
            D = InOutGenerator(_ffq,0.0)
            psi = build_psi(D)
            JSON.write(
                (@__DIR__)*"/data/psi/order_"*string(k1)*"_model_"*string(k2)*".json",
                JSON.json(psi),
            )

            # evaluate first-return distribution
            first_return_coeffs = zeros(size(models[k1][k2]["coeffs_0"]))
            first_return_coeffs[index(_ffq,Minus)] = 
                transpose(models[k1][k2]["coeffs_0"][index(_ffq,Plus)])*psi
            first_return_dist = SFMDistribution(first_return_coeffs,_ffq.dq)
            cdf_fun = cdf(first_return_dist)
            first_return_cdf_phase_2 = cdf_fun.(x_vals,2)
            first_return_cdf_phase_4 = cdf_fun.(x_vals,4)
            CSV.write(
                (@__DIR__)*"/data/first_return_cdf_approximations/order_"*string(k1)*"_model_"*string(k2)*".csv",
                DataFrame(
                    x = x_vals,
                    phase_2 = first_return_cdf_phase_2,
                    phase_4 = first_return_cdf_phase_4,
                ),
            )
        end
    end
end
