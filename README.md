### A barotropic model for flows in rotating-tank experiment
This is a finite-differencing model for the stream function-relative voticity equation:
$$\frac{\partial}{\partial t}(\nabla^2-k_d^2)\psi+J(\psi,\nabla^2\psi)+\beta\frac{\partial \psi}{\partial x}+\gamma\nabla^2\psi=Q$$
1. $\psi$ is the stream function whose dimension is $U^2L$ with $U$ and $L$ denoting the characteristic velocity scale and length scale;
2. $\zeta=\nabla^2\psi=(\partial_{xx}+\partial_{yy})\psi$ is the relative vorticity with a dimension of $1/T$ and time scale $T$;
3. $k_d=1/R_d$ is the reciprocal of the deformation radius $R_d=\sqrt{gH}/f_0$ with the gravitational acceleration $g$, mean water depth $H$, and the Coriolis parameter $f_0$;
4. $J(\psi,\zeta)=\psi_x\zeta_y-\zeta_x\psi_y$ is the nonlinear advection written in a Jacobian form;
5. $\beta$ is the topographic beta effect, it has a dimension of $1/(LT)$;
6. $\gamma$ is the Ekman decaying rate with dimension $T^{-1}$;
7. $Q$ is the forcing in terms of $\zeta$ injection thus its dimension is $T^{-2}$. 

The $\psi-\zeta$ equation has been historically employed to simplify the [Shallow Water Equations](https://en.wikipedia.org/wiki/Shallow_water_equations). In practice, an additional Laplacian visocisity $A_h \nabla^2\zeta$ is needed in RHS to damp the short waves generated by the nonlinear term and to maintain the numerical stability, where $A_h$ is the kenimatic vicosity with dimension $UL$. 

In this project, we use the $\psi-\zeta$ equation to model the barotropic component of the flow in the rotating tank experiments. The tank has a diameter of 110 cm and a depth of 40 cm; it is installed in a turnable table. Before the experiment, a solid barrier was placed inside the tank and was aligned in the raidal direction. This barrier acts as the west wall to the flows which was generated by a localized source placed to the East(azimuthal direction). 

The geometry configuration was captured by a camera as real-color pictures saved in `./for_grid_gen/`. The location/orientation of the wall was varied in some experiments, so you may find multiple images there. The following is an example showing one configuration in the experiment:

<img src="https://github.com/zy-gdou/rotating-tank-experiment/blob/385a69e7659770761996ae03dda59c8295bc769b/for_grid_gen/4oclock.png" width=50% height=50%>

Each image can be used by `genXYgrid_vect.m` to create a gird file for a specific geometry configuration. The grid file is named as `N?.mat` with `?` standing for the horizontal grid resolution(default value 500. It should be noted that if changing the resolution to a higher value, say 1000, the number of iteration used for the Possion-equation solver needs to be increased as well, otherwise the level of numerical convergence is low and the large-scale pattern cannot be presented correctly. Ps the computation cost for the high-resolution runs could be much more expensive than those low-resolution runs). The grid file must be generated before runing the code and is saved in the same directory as `genXYgrid_vect.m`. 

The square domain occupied by the tank is discretized into a horizontal Cartesian grid. Only the region inside the tank is modeled.
`genXYgrid_vect.m` defines the boundary points(the points located at the side wall of the tank), ghost points(points dropping outside of the tank with  a distance of 1 $dx$ to the neighouring boundary points, where $dx$ is the grid size), and inside points(points dropping inside the tank accounting for the major part of the total node population). This classification is useful when implementing the boundary conditions. Two boundary conditions: no-slip and free-slip boundary condtions can be applied to the lateral boundaries using the schemes in the review paper by E, Weinan & Liu, Guoqiang 1996(Journal of Computational Physics). The lateral boundary includes the side wall of the cylindrical tank and 3 sides of the barrier; the barrier is assumed to be a thin rectangle in `genXYgrid_vect.m` with one side joining the tank wall. Therefore there are 4 sides for the numerical domain. 

`sweep_run.m` defines the parameter space to survey and runs the model using each combination of these parameters. The key parameters are:
1. amplitude of the forcign vortex: Amp, $s^-2$ ;
2. radius of the forcing vortex: forcing_R , cm; 
3. barotropic deformation radius: Rd,cm; 
3. Ekman friction: gamma, $s^{-1}$; 
4. bulck viscosity: Ah, $cm^2/s$ 

Given the simplicity of the $\psi-\zeta$ equation, its dimensional form is discretized; this makes the code more readable and results easier to interpret. In `sweep_run.m`, user can also specify
1. the total time of integration:T_tot,s;
2. output frequency:dt,s; 
3. logical switch for the types of the boundary condition: noslip=true by default;
4. logical switch for nonlinear advection: nonlinear=true by default; otherwise the nonlinear term won't be included.

`pars.m` checks if this is a new run or a restart from any previous run. If there is no folder with the name specified as the parameter combination given in `sweep_run.m`, then a new folder will be created and the output files from the new run will be saved there. If there exists a folder with the same name as specified by the parameter combination in `sweep_run.m` and that folder contains multiple output files (.mat), then this is a restart run. Its initial condition is given by the last output file from the previous run. The restart run then append output files to the old folder. `pars.m` also contains other experimental parameters needed for the simulation and useful coefficients for the Possion equation solver and those for the time marching schemes. 

`main_iter_4step_RungeKutta.m` is the main iteration. It marches the governning eq. in time using 4th-order explicit Runge-Kutta scheme for the nonlinear term(with Arakawa 1966's conserving scheme for the spatial discretization) if nonlinear=true, and explicit (2nd-order central differencing) for the other terms(beta term, Ekman term, diffusion term). 

`allocate_matrices.m` initialize the useful matrices by allocating 1D vectors for them. This allocation is called everytime when the model starts (no matter if it is a new run or a restart run). When time integration finishes without overflow, the final-step $\psi$ will be saved in the case folder, which is a checkpoint for restart.

One example of the simulation results is shown in the [![video]()](https://youtu.be/MHM2IbPPq_U).

# Double-Fourier theory 
`bplume_slantwall_unrotate_lab.m` formulated a Double-Fourier(linear) theory. This theory is used to decompose the total wave field into the incidental and the reflected ones in a rectangular domain (a Cartesian coordinate system). It shows that the formation of the meanders at the flanks of the $\beta$-plume can be caused by the reflected Rossby waves.

