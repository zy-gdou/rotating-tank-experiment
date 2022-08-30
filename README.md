# rotating-tank-experiment
This is a finite-differencing model for the stream function-relative voticity equation used to model the barotropic flow in the rotating tank experiment 
performed at Memorial University of Newfoundland in 2022.

In order to use this model, a gird file ('N?.mat', where ? stands for the horizontal grid resolution) needs to be generated using "genXYgrid_vect.m". It creates a horizontal Cartesian grid covering the tank geometry(110cm in diameter). During the lab experiment, a solid wall was inserted into the tank and was aligned in the radial direction. This geometry configuration was recorded by a camera as full-color pictures saved in saved in 
'./for_grid_gen/' 

The location of the wall was varied in some experiments, so you may find multiple images there. Each image can be used by "genXYgrid_vect.m" to create the Cartesian grid with different resolutions(default resolution is 500)."genXYgrid_vect.m" also defines the boundary points and the ghost points, and insider points used to implement the boundary conditions(no-slip and free-slip) according to E Weinan & Liu,Guoqiang 1996.


'sweep_run.m' is the main code. User can speicify the forcing amplitude(Amp), radius of the forcing vortex(foring_R), deformation radius(Rd), 
Ekman friction gamma, viscosity Ah, total time of integration T_tot, and output frequency dt.

'pars.m' set the experimental parameters used in the simulation.

'allocate_matrices.m' initialize the matrices used by the equation; these 1D vectors are initialized everytime (regardless of restart or not) 
before the model runs.

'main_iter_4step_RungeKutta.m' marches the governning eq. in time using 4th-order explicit Runge-Kutta scheme for the nonlinear term(with Arakawa 1966's conserving 
scheme for the nonliner Jacobian), and implicit (and 2nd-order central differencing) for the other terms(beta term, Ekman term, diffusion term). 

