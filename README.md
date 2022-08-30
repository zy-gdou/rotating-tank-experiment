# 1.5-layer model for rotating-tank experiment
This is a finite-differencing model for the stream function-relative voticity equation used to model the barotropic flow in the rotating tank experiment with a cylindrical tank of 110cm diameter 110cm and 40cm depth . Before the experiment, a solid barrier inserted inside acting as the wall and it is aligned in the raidal direction.

In order to use this code, a gird file ('N?.mat', where ? stands for the horizontal grid resolution) needs to be generated using "genXYgrid_vect.m". It creates a horizontal Cartesian grid covering the tank geometry. This geometry configuration was recorded by a camera as full-color pictures saved in saved in './for_grid_gen/' 

The location of the wall was varied in some experiments, so you may find multiple images there. Each image can be used by "genXYgrid_vect.m" to create gird file with different resolutions(default horizontal resolution is 500)."genXYgrid_vect.m" also defines the boundary points, the ghost points, and insider points needed to implement the boundary conditions(no-slip and free-slip) according to E Weinan & Liu,Guoqiang 1996.

'sweep_run.m' is the main code. User can speicify the forcing amplitude(Amp), radius of the forcing vortex(foring_R), deformation radius(Rd), 
Ekman friction (gamma), viscosity (Ah), total time of integration (T_tot), and output frequency (dt) there.

'pars.m' contains other experimental parameters used in the simulation, and coefficients for the Possion equation solver and time marching.

'allocate_matrices.m' initialize the matrices used by the equation; they are 1D vectors initialized everytime (regardless of restart or not) 
before the model runs.

'main_iter_4step_RungeKutta.m' is the main iteratoin part. It marches the governning eq. in time using 4th-order explicit Runge-Kutta scheme for the nonlinear term(with Arakawa 1966's conserving scheme for the nonliner Jacobian), and implicit (and 2nd-order central differencing) for the other terms(beta term, Ekman term, diffusion term). 

