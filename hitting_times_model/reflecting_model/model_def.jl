# include("../preamble.jl")
include((@__DIR__)*"/../../preamble.jl")

T_reflecting = [-1.1 1.1; 1.0 -1.0]
c_reflecting = [1.0; -1.0]
b_reflecting = 10.0

P_upr_reflecting = [0.0 1.0]
P_lwr_reflecting = [1.0 0.0] # both reflecting boundaries

model = BoundedFluidQueue(T_reflecting,c_reflecting,P_lwr_reflecting,P_upr_reflecting,b_reflecting)

function build_discretised_model(type::Type{<:Mesh{T}},model,hx,order) where T 
    nodes = 0.0:hx:10.0
    mesh = type(nodes,order)
    dq = DiscretisedFluidQueue(model,mesh)
    B = build_full_generator(dq)
    d0 = interior_point_mass(eps(),1,dq)
    return d0, B, dq
end
t = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
dg_model(order) = build_discretised_model(DGMesh{t},model,1.0,order)
order1_model(order) = build_discretised_model(DGMesh{t},model,1.0/order,1)
qbdrap_model(order) = build_discretised_model(FRAPMesh{t},model,1.0,order)

