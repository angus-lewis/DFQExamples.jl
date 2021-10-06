**Function reconstrcution**
+ when DG of FV work, they work well (polynomial)
+ when they fail, they fail in bad ways (step function)
+ QBD-RAP always works but when the solution is smooth, pdf reconstruction is wiggly (polynomial)
+ all methods can reconstruct cell probabilities exactly

**One D wave** 
Scenario: start with point mass at 0, which moves to the right at rate 1, evlove to t for t = 1.5+eps(), 1.75, 2.0-eps(), (i.e. the left hand side, middle, and right hand side of a cell)
+ reconstruction of CDF/PDF
    + when we evolve over time, the wiggliness of DG and FV are accentuated, QBD-RAP works
+ cell averages
    + to remove the effects of reconstruction, we look at the probability of being in a cell. DG and FV have negative probs, sometimes drastically depending on where the true point mass lies within the cell.

**Variables**
+ Initial condition
+ Fluid queue
    + generator
    + rates
    + boundary conditions
+ Performance measures
    + transient distribution
    + stationary distribution
    + hitting times
    + first return of SFFM
+ space integration
    + methods
        + DG
            + first order, positivity preserving with cell width as parameter
            + high order
                + vanilla
                + flux limiters (positivity preserving)
                + post-hoc filtering (positivity preserving)
        + FV
            + high-order polynomial interpolation
            + filtering/limiting? (not sure what exists in this literature)
            + same approximation size? (not sure what exists in this literature)
        + QBD-RAP
            + closing vectors x3
            + could also use polynomial reconstructions on cell averages
    + cell sizes
+ time integration
    + scheme
        + RK4
        + Euler
    + step size
+ error metrics
    + cell probabilities (no reconstruction required)
    + performance of reconstruction (L^2 metric?)
+ second order models and pdes, MMBM

**Numerics chapter working outline**
+ Function reconstruction and closing operators
    + investigate order
    + functions:
        + point mass 
            + left cell edge
            + right cell edge
            + cell center
        + piecewise constant
            + discontinuities at:
                + left cell edge
                + right cell edge
                + cell center
        + piecewise linear
            + discontinuities at:
                + left cell edge
                + right cell edge
                + cell center
        + global quadratic
        + global trig/exponential
        + global ME
    + methods:
        + DG/polynomial projection with GLJ nodes
        + FV/polynomial projection with evenly spaced nodes
        + ME:
            + unnormalised 
            + niave normalisation 
            + geometric normalisation
+ Transient distributions and time integration
    + investigate order x t-step size
    + models:
        + 2 x fluid queues 
        + 1 x travelling wave
    + time integration: 
        + RK4
        + Euler
        + t-step size:
            + cell_width/max(c_i) * (0.1:0.1:1.0)
    + initial conditions:
        + point mass:
            + t small, corresponding to cell edge
            + t small, corresponding to cell centre
                + plot error vs t-step size for each order on same axis, arrange as 2d array of plots for time and space integration method, reproduce for each model x t
        + polynomial/exponential
            + t medium
                + plot error vs t-step size for each order x time integration method on same axis, arrange as 2d array of plots for time and space integration method, reproduce for each model
        + point mass at a boundary 
    

transient distributions loop variables 
    model: 2 x fluid queue + 1 x 1d wave                3
    integration method: RK4, Euler                      2
    space integration: DG, FV, QBD-RAP                  3
    initial condition: point mass, smooth, boundary     3
    t: small (mid-cell), small (cell edge), big         3
    cell sizes: n_cells = 1,2,4,8,16,24,48              7
    order: 1,3,5,7,9,15,31                              7
                                                        7938
save coeffs, generator

transient distributions simulation (n_sims = 10_000_000)
    model: 2 x fluid queue                              2
    initial condition: point mass, smooth, boundary     3
    t: small (mid-cell), small (cell edge), big         3
                                                        18
save all X(t), phi(t) values
