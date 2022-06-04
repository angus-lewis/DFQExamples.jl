# include((@__DIR__)*"/../preamble.jl")
function make_model() 
# generic parameters
γ₁ = 11.0; β₁ = 1.0; λ₁ = 12.48; θ₁ = 1.6; κ = 2.6
γ₂ = 22.0; β₂ = 1.0; λ₂ = 16.25; θ₂ = 1.0; xˢ = 1.6

# generator 
T = [
    -(γ₁+γ₂)        γ₂          γ₁         0.0
        β₂       -(γ₁+β₂)      0.0          γ₁
        β₁          0.0     -(γ₂+β₁)        γ₂
        0.0         β₁          β₂       -(β₁+β₂)
]

# first fluid rates, cᵢ
C = [
    λ₁-θ₁
    λ₁-θ₁
    -θ₁
    -θ₁
]

# second fluid rates, rᵢ(x)
# lower boundary
r_lwr = [
    # λ₂-κ
    # -κ
    λ₂-κ
    -κ
]

# in the interval (0,xˢ)
_size1(h) = (1,Int(round(1.6/h)))
r_below_xˢ(h) = [
    fill(λ₂-θ₂,_size1(h)...)
    fill(-θ₂,_size1(h)...)
    fill(λ₂-θ₂,_size1(h)...)
    fill(-θ₂,_size1(h)...)
]

# in the interval [xˢ,48) (truncate the model at 48)
_size2(h) = (1,Int(round(48.0/h-1.6/h)))
r_above_xˢ(h) = [
    fill(λ₂,_size2(h)...)
    fill(0.0,_size2(h)...)
    fill(λ₂,_size2(h)...)
    fill(0.0,_size2(h)...)
]

r_interior(h) = [r_below_xˢ(h) r_above_xˢ(h)]

# rates at the upper boundary, 48, due to truncation
# choose some sensible values: same as rates just below boundary
r_upr = [λ₂; 0.0] 

b = 48.0
first_fluid = BoundedFluidQueue(T,C,b)
model = first_fluid

h = 0.4
nodes = 0.0:0.4:b

ffq(mesh,h) = FluidFluidQueue(DiscretisedFluidQueue(first_fluid,mesh),r_interior(h),r_lwr,r_upr)

# simulation can be simplified if we use a coarser discretiseation
r_interior_sim = [r_below_xˢ(0.4)[:,1] r_above_xˢ(0.4)[:,end]]
# the true model is unbounded, so choose the upper bound of 
# the model we simulate from to be large enough that it 
# shouldn't affect the simulation dynamics
big_b = 1_000.0 
nodes_sim = [0.0;xˢ;big_b]

first_fluid_sim = BoundedFluidQueue(T,C,big_b)
mesh_sim = DGMesh(nodes_sim,1)
dq_sim = DiscretisedFluidQueue(first_fluid_sim,mesh_sim)

ffq_sim = FluidFluidQueue(dq_sim,r_interior_sim,r_lwr,r_upr)
return ffq_sim
end

const model = make_model()
