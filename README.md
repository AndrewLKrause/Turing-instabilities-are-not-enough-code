# Turing Conditions are not Enough Code
 
This is the repository of code associated with the paper, "Turing conditions are not enough to ensure pattern formation." The key files to use are RunAndPlot and the associated plotting functions. An example of this would be calling:

`[U, x, T, ui, vi] = RunAndPlot('KellerSegel',1);`

which will run the Keller-Segel chemotaxis model using the standard parameters in one spatial dimension, and then plot a Kymograph of the variable $u$.

All other functions are used to generate the Figures and Tables in the paper, such as the SystematicRuns.m file which implements the large-scale Latin Hypercube Sampling and runs the corresponding simulations.

## Numerical methods

For the three local in space models (reaction-diffusion, Keller-Segel, and the Biharmonic model), a simple central finite differencing method is used to discretize the spatial operators. For the nonlocal advection model, we only include a 1D code here given in the NonlocalAdvectionSolver.m which implements a pseudospectral discretization of the integral. 2D and 1D codes corresponding to that model can be found at https://github.com/JunJewell/NonlocalReactAdvectDiffuse2D with further documentation.

In all 1D simulations we used $m=1000$ nodes, and for the 2D simulations we used $m=100^2$ nodes (with $m=300^2$ nodes used to produce Figure 3). In all of these cases, the spatially-discretized system is then fed into the MATLAB function `ode15s` which implements a stiff implicit timestepping scheme, with default relative and absolute tolerances set to $10^{-9}$. For the local models, Jacobian sparsity patterns are used to massively reduce simulation times, though this is not useful for the nonlocal model which generically has a fully-coupled system of ODEs after spatial discretization. Convergence checks were carried out in space and time for the standard parameters of each model.
