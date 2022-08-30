% linear beta-plume reflected at a slanted wall placed at the west end (Xwall=0)
% modelled using a Double-Fourier decomposition.

% The script gives the solutions for the incidental plume, reflected plume, and their
% sum in terms of surface elevation, at a series of moments,
% the parameters used here are for the rotating-tank exp.

% The coordiante sys is the original x and y axes, with x pointing to the east and y to the
% north. The wall starts from x=0,y=0, and slanted with an angle, theta, with respect to y
% pointing to the northeast,
clear all;
close all;
clc

% nv=3;   % for quiver plot
% uv_scale=3;

moments=[0.1:0.1:70] ;   % in second
% ep =tan(2*pi*1.5/12*0.5);           % tan(alpha), alpha is the angle between the wall and y axis,
ep =tan(pi*26.6/180);           % tan(alpha), alpha is the angle between the wall and y axis,
alpha_angle =  atan(ep); % angle between the wall and y axis

Amp = 0.1;         % s^-2, magnitude of Q which foces the vort  eq. due to the eddy forcing
forcing_R = 4.5;    % cm, radius of the Gaussian vortex formed by the injected freshwater;
gamma = 0.025;    % sqrt(0.041*2.6)/8;     % 1/(8e1); % <= Ekman friction coefficient= sqrt(\nu*OM)/H
beta=0.12; % 2*OM^3*35/(g*H); %topographic beta effect, 1/(s*cm)
Rd =4.5; %forcing_R;
kds = Rd^(-2) ;%; %<=kd square =Rd^-2, Rd is the b.c. deformation radius


Ly=55*2; %cm; radis of the tank or the barrier length
Lx=2*pi*35*5/12;% pi*Ly/2; %zonal extent of the beta-plume obseved in the exp. 2*pi*35*5/12*2
Nx=128; % no. of grid point in the zonal direction
Ny=128; % no of grid point in the south-north direction
forcing_x=0.72;%ep*0.5+0.5/cos(theta); % Nx*0.5; %distance to the western boundary
forcing_y=0.5; % distance to the souther/northern boundary

g=982; %cm/s^2
OM=2.3; %2*pi/sec_ped; % rad/s
f=2*OM;%*sind(lat);    % rad/s
H=10;                  % mean water depth, cm

% gprime=g*1.e-3; %reduced gravity
% H1=H*1.e-1;            % upper layer thickness if necessary


% define the horizonatal grid
[X,Y]=meshgrid(0:Nx-1,0:Ny-1); %horizontal grid. unit:dx=Lx/Nx, dy=Ly/Ny
X=X/Nx*Lx;
Y=Y/Ny*Ly;
dx=Lx/Nx;dy=Ly/Ny;
Xwall=0;
%X(1,Nx*forcing_x+30); %location of the wall, must to the east of the forcing

yb=(0:Ny-1)/Ny*Ly; %cooridinates of the wall
xb=ep*yb;





% define the Gaussian forcing Q, which drives the vort equation, its dimension is s^-2
Q=Amp*exp(  -(X-forcing_x*Lx).^2./forcing_R^2 ....
    -(Y-forcing_y*Ly).^2./forcing_R^2 ); %dimension, 1/s^2
fQ=fft2(Q)/(Nx*Ny);   % which enters into the solution

%  define the wavenumber space
kx=[0:Nx/2,-Nx/2+1:-1]*2*pi/Lx; %cm^-1
ky=[0:Ny/2,-Ny/2+1:-1]*2*pi/Ly;
[Kx,Ky] = meshgrid(kx,ky);
ks  = Kx.^2+Ky.^2;
if gamma*kds~=0
    Kx_r=beta*(Kx.^2+Ky.^2+kds)*cos(alpha_angle)^2;
    Kx_r=Kx_r./(Kx*beta-1j*gamma*kds);
    Kx_r=Kx_r+Ky*sin(2*alpha_angle)-Kx*cos(2*alpha_angle);
else
    Kx_r= (Ky.^2+kds)./Kx ;
end

%dispersion relation of the Rossby wave, eq.(8), s^-1
omega= ( -Kx.*beta -1j*(gamma*ks) )./(ks+kds);
omega(isnan(omega))=0;

Ky_r=ep*(Kx-Kx_r)+Ky;


% since d\zeta/dt = Q*g/f, Q*g/f*dt gives the relative vorticyt increment
% in dt time injected by the forcing, here dt=1day=86400s, then Q*g/f^2*86400s gives the
% Ro number increment after being forced by 1day. Therefore, this figure
% quantifies how strong the forcing is
% figure;
% pcolor(X,Y, Q*g/f^2);
% xlabel('x,km');ylabel('y,km')
% title( 'Ro number increment after being forced by 1s' );
% colorbar;hold on
% % plot(xb*1.e-3,yb*1.e-3,'w-','linewidth',5);
% text(20 ,400 ,[num2str( 90-atan(1/ep)/pi*180,'%3.0f'),'^o'],'fontsize',15,'color','w')
foldername=['Amp=',num2str(Amp,'%3.1f'),' Rf=',num2str(round(forcing_R)),.....
    ' Rd=',num2str(round(Rd)),' frc=',num2str(gamma,'%4.3f'),....
    ' alpha=',num2str(alpha_angle/pi*180,'%3.0f')];
[status,msg]=mkdir(foldername);
if ~isempty(msg)
    disp([foldername,' existed!'])
else
    disp([foldername,' is created! pictures are saved there.'])
end
movie_name=['theory movie ',num2str(moments(1)),'-', num2str(moments(end)),'s'];
    
% v=VideoWriter([foldername,'/',movie_name]);
% v.FrameRate=50;
% v.Quality=100;
% open(v);


% fl=dir([foldername,'/*.png']);
% [~,idx]=sort([fl.datenum]);
% fl=fl(idx);
% for it=1:numel(fl)
% im=imread([foldername,'/',fl(it).name]);
% [imind,cm]=rgb2ind(im,256);
% if it==1
% imwrite(imind,cm,[foldername,'/gifmovie.gif'],'gif','DelayTime',0.1,'Loopcount',inf);
% else
%     imwrite(imind,cm,[foldername,'/gifmovie.gif'],'gif','DelayTime',0.2,'WriteMode','append');,
% end
% 
% end
  
wall=X<xb'*ones(1,Nx);
wall2=X==xb'*ones(1,Nx);



%drawnow;
% for it=1:numel(moments)
%     t=moments(it);
t=36;
    T_rot=num2str(t/(2*pi/OM),'%3.0f');
         disp([' @ ',num2str(t), 's, ',T_rot,' rot.']);
         

                  
    % eq.(12) in the doc file gives the incident plume without the boundary effect
    ss= fQ./(-1j*beta.*Kx+gamma.*ks).*(1-exp(-1j*omega*t))*Nx*Ny;
    ss(isnan(ss))=0;
    psi_inc = real(ifft2(ss));
%     u=real(ifft2(ss.*1j.*-Ky)); %m/s
%     v=real(ifft2(ss.*1j.*Kx));  %m/s
    
    % the total psi as  sum of incident and reflect
    ss=zeros(size(X));    
    for n2=1:Nx*Ny
        deno =-1j*beta*Kx(n2)+gamma*ks(n2);
        if ( deno~=0 && ~isnan(Kx_r(n2)) && ~isinf(Kx_r(n2)) )
            T = fQ(n2)/deno*(1-exp(-1j*omega(n2)*t));
            ss=ss + T.*(.......
                exp(1j*( Ky(n2)*Y   + Kx(n2)*X   )).......
                -exp(1j*( Ky_r(n2)*Y + Kx_r(n2)*X )).......
                );            
        end
    end
    psi=real(ss);
    psi(wall)=nan; %assign nan to the continent
    psi(wall2)=0;
    
%      zeta= ( psi(1:end-2,2:end-1)+psi(3:end,2:end-1) - 2*psi(2:end-1,2:end-1) )/dy^2+....
%            ( psi(2:end-1,1:end-2)+psi(2:end-1,3:end) - 2*psi(2:end-1,2:end-1))/dx^2 ;
    
    psi_ref=psi-psi_inc;
%    psi_inc(wall)=nan; %assign nan to the continent
    
    
    
    
    
    
    
    
    
    
    
    figure('position',[100,0,600,920],'color','w','Inverthardcopy','off');

% ========== draw these beta plumes and their sum ==============

    
   h1=axes('position',[0.12,0.7-0.02,0.82,0.28]);
    contourf(X,Y-Ly/2,psi*f/g,15);
   set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 75],'ylim',[-1 1]*30,'fontsize',12)
    hb=colorbar;
    hb.Color='k';
    hb.FontSize=12;
    caxis([-0.005 0.03])
    hold on
% text(2,25,[ T_rot,' rot.'],'color','w','fontsize',18,'fontname','serif');
% text(75,25,['\eta'],'color','w','fontsize',25,'fontname','serif');
% ylabel('y, cm','fontsize',18)
 

slice=psi_inc*f/g;
slice(wall)=nan;
   h2=axes('position',[0.12,0.7-0.28-0.06,0.82,0.28]);
    contourf(X,Y-Ly/2,slice,8);
        set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 75],'ylim',[-1 1]*30,'fontsize',12)  
     hb=colorbar;
    hb.Color='k';
    hb.FontSize=12;
%     ylabel('y, cm','fontsize',18)
    caxis([-0.005 0.03])



   h3=axes('position',[0.12,0.69-0.28*2-0.09,0.82,0.28]);
    contourf(X,Y-Ly/2,psi_ref*f/g,8);
    hb=colorbar;
    hb.Color='k';
          caxis([-6 2]*1.e-3)

        set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 75],'ylim',[-1 1]*30,'fontsize',12)
%  text(68,25,['\eta_{ref}'],'color','w','fontsize',25,'fontname','serif');
              hb.FontSize=12;
%     xlabel('x, cm','fontsize',18);
%     ylabel('y, cm','fontsize',18)

    print(gcf,[foldername,'/t=',num2str(t),'s.png'],'-r200','-dpng')
    
    
    
    

    
    
    
    
    
        
    figure('position',[100,0,300,920],'color','w','Inverthardcopy','off');    
   h1=axes('position',[0.13,0.7-0.02,0.82,0.28]);
    contourf(X,Y-Ly/2,psi*f/g,20);
   set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 8],'ylim',[-1 1]*30,'fontsize',12)
    caxis([-0.005 0.01])
    hold on
 colorbar;
 
   h2=axes('position',[0.13,0.7-0.28-0.06,0.82,0.28]);
    contourf(X,Y-Ly/2,slice,10);
        set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 8],'ylim',[-1 1]*30,'fontsize',12)  
    caxis([-0.005 0.005])
 colorbar;

   h3=axes('position',[0.13,0.69-0.28*2-0.09,0.82,0.28]);
    contourf(X,Y-Ly/2,psi_ref*f/g,8);
          caxis([-0.01 0.005])
          colorbar;
        set(gca,'ydir','normal','color','k','xcolor','k','ycolor','k',......
        'xlim',[0 8],'ylim',[-1 1]*30,'fontsize',12)
    print(gcf,[foldername,'/t=',num2str(t),'s near wall.png'],'-r200','-dpng')
    
    
    
    
    


