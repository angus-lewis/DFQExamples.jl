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
