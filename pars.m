ReS=501;
% note that when ReS doubles, the no. of iteration need in the SOR solver
% needs to be increased as well, which slows the simulation

% load the grid file
load (['N',num2str(ReS),' BarrierOiren 4.0.mat']);

x=linspace(0,Lx,N);  %x-grid used for drawing pictures
[Y,X]=meshgrid(x,x); %horizonal grid
r=sqrt((X-Lx*0.5).^2+(Y-Lx*0.5).^2); %radial distance to the axis of rotation, it is used to identify the boundary of the domain
dx=Lx/(N-1);   %spatial resolution, cm
rdx=1./dx;
rdxs=rdx*rdx;
dxs=dx^2;


% horizontal diffusion is needed only for the nonlinear run in orde to
% filter out the small waves.
% Amp= [-3]*1.e-1 ; %s^-1,
% Ah = [0.06]; %water kinematic viscosity at 20C deg. cm^2/s
% Ah=5e-1; %cm^2/s, horizonal viscosity for the Laplacian diffusion


OM=2.3;           % rad/s
g=980;           % cm/s^2, gravitational acceleration
A=OM^2/(2*982);  %controls the shape of the paraboloidal surface
f_Coriolis=2*OM; % Coriolis parameter in our exp
H0=10;           % cm, mean water depth
% beta=8.e-2; %1/(s*cm), topographic beta effect used in our lab experiment
Urms=1.5 ;%cm/s % typical max horizontal velocity observed in the experiment
dt= max([min([ dxs/Ah*0.5, 2*Ah/Urms^2])*0.5 0.02]);  %time step, s, E weinan and Liu guoqiang 1996 jcf
Nnt=round(outfreq_dt/dt); %every Nnt*dt seconds draw the flow field;

% T_tot=120;
rdxs_Ah=rdx*rdx*Ah;
dtby6=dt/6;
dtby2=dt*0.5;

H = H0+A*(r(ij(1:Ninsider)).^2-0.5*(Lx*0.5)^2);%parabolic surface
beta_T = f_Coriolis*2*A*r(ij(1:Ninsider))./H; %topographic beta effect
coscita = (X(ij(1:Ninsider))-0.5*Lx)./r(ij(1:Ninsider));
sincita = (Y(ij(1:Ninsider))-0.5*Lx)./r(ij(1:Ninsider));
coscita(isnan(coscita))=0; %in case that the polar point gives nan value
sincita(isnan(sincita))=0;

% Ekman friction term
% Ek_damp = 1/(1+8e-3*dt);%Ek friction coeff is 1e-2
%coefficiant gamma for the Ekman damping term in the relative vort eq.
% Ek_gamma=8.e-3; %s^-1
% Ek_gamma = sqrt(OM*Ah)./H; %shallow water decays faster than deep-water region
% Ek_gamma = median(Ek_gamma)*10;
Ek_lateral = 0.5; %lateral Ekman friction for the sidewall and the west barrier when non-slipe b.c. implemented

% forcing_R = 6; %cm, radius of the Gaussian forcing
sd2sw=35; %cm, radial distance from the source to center

%External Guassian forcing act on the RHS of the vorticity eq
%The following forcing term goes into the RHS of the vorticity equation.
Qforcing = Amp*exp(.....
    -((X(ij(1:Ninsider))-(Lx*0.5+real(nv2source)*sd2sw)).^2+(Y(ij(1:Ninsider))-(Lx*0.5+imag(nv2source)*sd2sw)).^2).....
    /forcing_R^2 .....
    );
% Qforcing's dimension is s^-2 ,[] means taking dimension of the quantity,
% then one has from eq.(6) of the draft that
% [g/f_Coriols][Freshwater_Volflux/(pi*forcing_R^2)][kds]=[Qforcing]=[Amp]=s^-2

% Rd=19;%sqrt(g*H0)/f_Coriolis; %cm, barotropic deformation radius
% kds=forcing_R^-2*kds_n;%  kd^2, kd=1/Rd,  is deformation radius
kds = Rd^(-2);% barotropic case

% estimated freshwater flux through the tube, cm^3/s
Freshwater_Volflux=abs(Amp*f_Coriolis/(g*kds) *forcing_R^2*pi) ; %max freshwater vol flux
%[Amp*f_Coriolis/g/kds]=cm/s,
% g/f_Coriolis*kds*   Freshwater_Volflux/(forcing_R^2*pi); % zeta increment
% in 1 second

lambda_SOR = 1.75; %Coefficiant for the Successive Over-relaxation method
Coe_SOR = lambda_SOR./(4+kds*dxs);
Coe_Jac = 1./12.0000000000000000000*rdx^2;% Coefficient used before the Arakawa's 9-point Jacobe operator
% Ro_inc = num2str(abs(Amp),'%4.1e'); %% Qforcing *dt can be used to quantify the strenth of nonlinearity
% clear X Y;

%  bottom line of Weinan E&Liu 1996's tabel 1: for no-slipe b.c.
OrszagIsraeli74_35o13= 35/13;
OrszagIsraeli74_m1o13= -1/13;
OrszagIsraeli74_22o13= 22/13;



% linear =false; % if ture, nonlinear Jacobe is turned off, otherwise nonlinear term is included
% no_slip=true; %true for no-slip boundary condition on 4 solid sides, false for free-slip b.c


% create a folder to save the output images and mat files from each run
if linear
    if no_slip
        subtitle='linear non-slip';
    else
        subtitle='linear free-slip';
    end
else
    if no_slip
        subtitle='non-linear non-slip';
    else
        subtitle='non-linear free-slip';
    end
end

case_folder=['Amp=',num2str((Amp),'%5.3f'),' R_f=',num2str(forcing_R),......
    'cm Rd=',num2str( Rd,'%3.0f'),.......
    'cm Ah=',num2str(Ah,'%3.2f'),' frc=',num2str(median(Ek_gamma),'%6.4f')   ];


dir2cf=['N' num2str(ReS),'/Ah=',num2str(Ah),'/'];
[status, msg]=mkdir([dir2cf,case_folder]);
if ~isempty(msg)
    clear fls
    fls=  dir([dir2cf,case_folder,'/*.mat']);
    if numel(fls)>1
        [~,idx] = sort([fls.datenum]);
        fls = fls(idx);
        disp([dir2cf,'/',case_folder,'"\ existed, & contains many mat files; '])
        T=str2num(fls(end).name(3:end-5));
        Niter=0;
        restart=true;
        disp(['Restart run (',subtitle,') from ',num2str(round(T)),' s to ', num2str(round((T_tot+T))),' seconds'])
        disp(['Output frequencey is ',num2str(Nnt*dt),' seconds.'])
        T_tot=T_tot+T;
        
    else
        disp([dir2cf,case_folder,'"\ existed, but contains no mat files; '])
        disp(['New run(',subtitle,') to ',num2str(T_tot),' seconds.' ]);
        disp(['Output frequencey: ',num2str(Nnt*dt),' seconds.'])
        T = 0;
        Niter = 0;
        restart=false;
    end
else
    disp([dir2cf,case_folder,'"\ is created; ']);
    disp(['New run(',subtitle,') to ',num2str(T_tot), ' seconds; ']);
    disp(['Output frequency is ',num2str(Nnt*dt) ,' seconds.'])
    restart=false;
    T=0; %total integration time, s
    Niter=0; %total number of time integration.
end


disp('----------------- Forcing parameters ---------------')
disp(['Forcing eddy size=',num2str(forcing_R,'%3.1f'),......
    'cm, Rd=',num2str(Rd,'%3.1f'),'cm;',......
    'Amp=',num2str(Amp),' Fw flux=',num2str(Freshwater_Volflux,'%5.3f'),' cm^3/s'])
disp(['---------------- Damping parameters----------------'])
disp([' Ah=',num2str(Ah,'%4.3f'),' cm^2/s; gamma=',num2str(median(Ek_gamma),'%6.1e'),'s^-1, dt=',num2str(dt,'%10.5f'),'s.'])

