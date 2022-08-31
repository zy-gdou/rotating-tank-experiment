# 1.5-layer model for rotating-tank experiment
This is a finite-differencing model for the stream function-relative voticity equation which is used to model the barotropic flow in the rotating tank experiment. The tank has a diameter of 110cm and a depth of 40cm and is installed in a turnable table. Before the experiment, a solid barrier was placed inside the tank and was aligned in the raidal direction. This barrier acts as the west wall to the flows generated by a localized source located to the east. 

The geometry configuration was recorded by a camera as full-color pictures saved in './for_grid_gen/'. Since the location/orientation of the wall was varied in some experiments, you may find multiple images there. Each image can be used by "genXYgrid_vect.m" to create a gird file for a specific geometry configuration. The grid file is named as "N?.mat" with "?" standing for the horizontal grid resolution(default value 500).
![This is an image](https://github.com/zy-gdou/rotating-tank-experiment/blob/main/for_grid_gen/1oclock.tif?raw=true "geometry setup in one experiment")
A grid file must be generated before runing the code and it is saved in the same directory as "genXYgrid_vect.m". 

Our simulation employes a horizontal Cartesian grid to cover the entire area captured by the camera, and only the region inside the tank is modeled.
"genXYgrid_vect.m" defines the boundary points(the points located at the side wall of the tank), ghost points(points dropping outside of the tank with  a distance of 1 dx to the neighouring boundary points, where dx is the grid size), and inside points(points dropping inside the tank which are the major part of the total node population). This classification is useful to implement the boundary conditions. Two boundary conditoins (no-slip and free-slip) can be applied in the code according to the review paper by E, Weinan & Liu, Guoqiang 1996(Journal of Computational Physics).

'sweep_run.m' defines the parameter space we are going to survey. The key parameters are:
the amplitude of the forcign vortex: Amp, s^-2 ;
the radius of the forcing vortex: foring_R , cm; 
barotropic deformation radius: Rd,cm; 
the Ekman friction: gamma, s^-1; 
viscosity: Ah, cm^2/s 
Given the simplicity of the stream function-relative vorticity equation, its dimensional form is discretized; this makes the code more readable and results easier to interpret.
In "sweep_run.m", user can also specify
the total time of integration:T_tot,s;
output frequency:dt,s; 
logical switch for the types of the boundary condition: noslip=true by default;
logical switch for nonlinear advection: nonlinear=true by default; otherwise the nonlinear advection term won't be included.

'pars.m' checks if this is a new run or a restart from any previous runs. If there is no folder with the name specified as the parameter combination given from sweep_run.m, then a new empty folder will be created, a new run starts, and the output files will be saved in that new folder. If there is an exist old folder with the same name as specified by the "sweep_run.m" and that folder contains multiple output files (.mat) inside, then this run will be recognized as an restart run. Its initial condition is given by the last output file. The restart run then adds output files to the old folder.
'pars.m' also contains other experimental parameters needed for the simulation and useful coefficients for the Possion equation solver and the time marching schemes. 

'main_iter_4step_RungeKutta.m' is the main iteration part. It marches the governning eq. in time using 4th-order explicit Runge-Kutta scheme for the nonlinear term(with Arakawa 1966's conserving scheme for the spatial discretization) if nonlinear=true, and explicit (2nd-order central differencing) for the other terms(beta term, Ekman term, diffusion term). 

'allocate_matrices.m' initialize the matrices for the variables used by the simulation by allocating 1D vectors for them. This allocation is called everytime when the model starts (no matter if it is a new run or a restart run) .

'bplume_slantwall_unrotate_lab.m' is the Double-Fourier(linear) theory used to decompose the total wave field into the incidental and the reflected ones in a rectangular domain (a Cartesian coordinate system). This theoy is used to show the formation of the meanders at the flanks of the beta-plume caused by the reflected Rossby waves.


