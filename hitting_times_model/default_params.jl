orders = 1:2:21#21
n_err_evals = 5_001
n_sims = 1_000_000 # need lots of these for good error plots
n_boot = 200
qtiles = (0.05,0.95)
transient_time = 2.0
approx_types = (:dg,:dg,:order1,:qbdrap) 
