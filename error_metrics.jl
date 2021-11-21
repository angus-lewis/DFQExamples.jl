kolmogorov_smirnov(cdf_1_evaluated, pm_1, cdf_2_evaluated, pm_2, grid) = maximum(abs.(cdf_1_evaluated-cdf_2_evaluated))
function L1(f1_evaluated, pm_1, f2_evaluated, pm_2, grid) 
    d = f1_evaluated-f2_evaluated
    d = sum(abs.(d[1:end-1,:] + d[2:end,:])./2.0.*diff(grid))
    return log10(d)
end
function L2(f1_evaluated, pm_1, f2_evaluated, pm_2, grid) 
    d = f1_evaluated-f2_evaluated
    d = sum(((d[1:end-1,:] + d[2:end,:]).^2)./2.0.*diff(grid))
    return log10(d)
end
function L1_cell_probs(f1_evaluated, pm_1, f2_evaluated, pm_2, grid) 
    d = sum(abs.(f1_evaluated-f2_evaluated))
    d += sum(abs.(pm_1[:]-pm_2[:]))
    return log10(d)
end
