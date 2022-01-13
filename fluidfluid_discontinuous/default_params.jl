orders = 1:2:21#; orders = [3;7;9;11;13;15;21]
n_err_evals = 10_001 # x-vales to evaluate the error metrics with
n_sims = 5_000_000 # need lots of these for good error plots
n_boot = 1_000 # number of bootstrap samples to perform to estimate confidence intervals
qtiles = (0.05,0.95) # 90% confidence intervals
transient_time = 2.0 # time to evolve the solutions to
approx_types = (:dg,:dg,:order1,:qbdrap) 
