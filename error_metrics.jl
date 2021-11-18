kolmogorov_smirnov(F1::Function,F2::Function,grid) = 
    maximum(abs.(F1.(grid)-F2.(grid)))

function kolmogorov_smirnov(F1::Function,F2::Function,grid,_phases) 
    M = 0.0
    for i in _phases
        m = maximum(abs.(F1.(grid,i)-F2.(grid,i)))
        (m>M)&&(M=m) 
    end
    return M
end

function kolmogorov_smirnov(d::SFMDistribution,F2::Function; n_evals=2001)
    grid = range(d.dq.mesh.nodes[1],d.dq.mesh.nodes[end],length=n_evals)
    return kolmogorov_smirnov(d,F2,grid)
end
function kolmogorov_smirnov(d::SFMDistribution,F2::Function,grid)
    F1 = cdf(d)
    M = 0.0
    for i in 1:n_phases(d.dq)
        m = maximum(abs.(F1.(grid,i)-F2.(grid,i)))
        (m>M)&&(M=m) 
    end
    return M
end
kolmogorov_smirnov(d::SFMDistribution,F2::Function,grid,_phases) = kolmogorov_smirnov(d,F2,grid)
function kolmogorov_smirnov(d::SFMDistribution,F2::Array,sim_pm::Array,grid,phases)
    F1 = cdf(d)
    M = 0.0
    n₋_counter = 0
    n₊_counter = 0
    for i in 1:n_phases(d.dq)
        m = max(maximum(abs.(F1.(grid,i)-F2[:,i])))
        if DiscretisedFluidQueues._has_left_boundary(model.S,i)
            n₋_counter += 1
            (abs.(F1(0.0,i) - sim_pm[n₋_counter])>m)&&(m=abs.(F1(0.0,i) - sim_pm[n₋_counter]))
        end
        if DiscretisedFluidQueues._has_right_boundary(model.S,i)
            n₊_counter += 1
            (abs.(F1(0.0,i) - sim_pm[end-N₊(d.dq.model)+n₊_counter])>m)&&(m=abs.(F1(d.dq.model.b,i) - sim_pm[end-N₊(d.dq.model)+n₊_counter]))
        end
        (m>M)&&(M=m) 
    end
    return M
end
kolmogorov_smirnov(d::SFMDistribution,s::Simulation,grid) = kolmogorov_smirnov(d,cdf(s),grid)
kolmogorov_smirnov(d::SFMDistribution,s::Simulation,grid,_phases) = kolmogorov_smirnov(d,cdf(s),grid)
kolmogorov_smirnov(d::SFMDistribution,s::SFMDistribution) = kolmogorov_smirnov(d,cdf(s))


function Lp_cell_probs(d::SFMDistribution,F2::SFMDistribution; p=1)
    grid = (d.dq.mesh.nodes[1:end-1]+d.dq.mesh.nodes[2:end])./2.0
    F1 = cell_probs(d)
    M = 0.0
    for i in 1:n_phases(d.dq)
        M += sum(abs.(F1.(grid,i)-F2.(grid,i)).^p)
    end
    M += sum(abs.(d.coeffs[1:N₋(d.dq)] - F2.coeffs[1:N₋(d.dq)]))
    M += sum(abs.(d.coeffs[end-N₊(d.dq)+1:end] - F2.coeffs[end-N₊(d.dq)+1:end]))
    return M
end

function Lp_cell_probs(d::SFMDistribution,F2::Simulation; p=1)
    grid = (d.dq.mesh.nodes[1:end-1]+d.dq.mesh.nodes[2:end])./2.0
    F1 = cell_probs(d)
    Fun2 = cell_probs(F2,d.dq.mesh.nodes)
    M = 0.0
    n₋_counter = 0
    n₊_counter = 0
    for i in 1:n_phases(d.dq)
        M += sum(abs.(F1.(grid,i)-Fun2.(grid,i)).^p)
        if DiscretisedFluidQueues._has_left_boundary(model.S,i)
            n₋_counter += 1
            sim_pm = sum(F2.X[F2.φ.==i].==0.0)/length(F2.t)
            M += abs.(d.coeffs[n₋_counter] - sim_pm)^p
        end
        if DiscretisedFluidQueues._has_right_boundary(model.S,i)
            n₊_counter += 1
            sim_pm = sum(F2.X[F2.φ.==i].==d.dq.model.b)/length(F2.t)
            M += abs.(d.coeffs[end-N₊(d.dq)+n₊_counter] - sim_pm)^p
        end
    end
    return M
end

function Lp_cell_probs(d::SFMDistribution,F2::Array,sim_pm::Array,grid,_phases; p=1)
    # grid = (d.dq.mesh.nodes[1:end-1]+d.dq.mesh.nodes[2:end])./2.0
    F1 = cell_probs(d)
    M = 0.0
    n₋_counter = 0
    n₊_counter = 0
    for i in _phases
        M += sum(abs.(F1.(grid,i)-F2[:,i]).^p)
        if DiscretisedFluidQueues._has_left_boundary(model.S,i)
            n₋_counter += 1
            M += abs.(d.coeffs[n₋_counter] - sim_pm[n₋_counter])^p
        end
        if DiscretisedFluidQueues._has_right_boundary(model.S,i)
            n₊_counter += 1
            M += abs.(d.coeffs[end-N₊(d.dq)+n₊_counter] - sim_pm[end-N₊(d.dq)+n₊_counter])^p
        end
    end
    return M
end

L1_cell_probs(d::SFMDistribution,F2::Simulation) = Lp_cell_probs(d,F2,p=1)
L1_cell_probs(d::SFMDistribution,F2::Array,sim_pm::Array,grid,_phases) = Lp_cell_probs(d,F2,sim_pm,grid,_phases,p=1)
L1_cell_probs(d::SFMDistribution,F2::Simulation,grid,phases) = L1_cell_probs(d,F2)

# Lp_cell_probs(d::SFMDistribution,s::Simulation,p=1) = Lp_cell_probs(d,cell_probs(s),p)
# Lp_cell_probs(d::SFMDistribution,s::SFMDistribution,p=1) = Lp_cell_probs(d,cell_probs(s),p)

function Lp_pdf(d::SFMDistribution,f2::Function,n_evals=2001;p=1)
    grid = range(d.dq.mesh.nodes[1],d.dq.mesh.nodes[end],length=n_evals)
    f1 = pdf(d)
    M = 0.0
    for i in 1:n_phases(d.dq)
        fn_evals = abs.(f1.(grid,i)-f2.(grid,i)).^p
        m = sum((fn_evals[1:end-1]+fn_evals[2:end])./2.0 .* diff(grid))
        M += m
    end
    return M
end
Lp_pdf(d::SFMDistribution,s::SFMDistribution;p=1) = Lp_pdf(d,pdf(s),p)

function Lp(f1::Function, f2::Function, grid; p=1)
    fn_evals = abs.(f1.(grid)-f2.(grid)).^p
    m = sum((fn_evals[1:end-1]+fn_evals[2:end])./2.0 .* diff(grid))
    return m
end

function Lp(f1::Function, f2::Function, grid, _phases; p=1)
    m = 0.0
    for i in _phases
        fn_evals = abs.(f1.(grid,i)-f2.(grid,i)).^p
        m += sum((fn_evals[1:end-1]+fn_evals[2:end])./2.0 .* diff(grid))
    end
    return m
end
function Lp(f1::Function, f2::Array, grid, _phases; p=1)
    m = 0.0
    for i in _phases
        fn_evals = abs.(f1.(grid,i)-f2[:,i]).^p
        m += sum((fn_evals[1:end-1]+fn_evals[2:end])./2.0 .* diff(grid))
    end
    return m
end

L1(f1::Function, f2::Function, grid, _phases) = Lp(f1, f2, grid, _phases; p=1)
L2(f1::Function, f2::Function, grid, _phases) = Lp(f1, f2, grid, _phases; p=2)
L1(f1::Function, f2::Array, sim_pm::Array, grid, _phases) = Lp(f1, f2, grid, _phases; p=1)
L2(f1::Function, f2::Array, sim_pm::Array, grid, _phases) = Lp(f1, f2, grid, _phases; p=2)
