# rotating-tank-experiment
Finite differencing model for the stream function-relative voticity equation used to model the rotating tank experiment

This is matlab code for the 1.5-layer stream function-relative vorticity equation employed to model the "barotropic" flow in the rotating tank experiment 
performed at Memorial University of Newfoundland in 2022.

A horizontal Cartesian grid is generated using the tank geometry(110cm in diameter). The geometry of the domain was recorded by a camera as full-color pictures during the experiment. 
These pictures are saved in 
'./for_grid_gen/' 
These pictures are used to generate grids of different resolution(default horizontal resolution is 500).

'sweep_run.m' is the main code. User can speicify the forcing amplitude(Amp), radius of the forcing vortex(foring_R), deformation radius(Rd), 
Ekman friction gamma, viscosity Ah, total time of integration T_tot, and output frequency dt.

'pars.m' set the experimental parameters used in the simulation.

'allocate_matrices.m' initialize the matrices used by the equation; these 1D vectors are initialized everytime (regardless of restart or not) 
before the model runs.

'main_iter_4step_RungeKutta.m' marches the governning eq. in time using 4th-order explicit Runge-Kutta scheme for the nonlinear term(with Arakawa 1966's conserving 
scheme for the nonliner Jacobian), and implicit (and 2nd-order central differencing) for the other terms(beta term, Ekman term, diffusion term). 

