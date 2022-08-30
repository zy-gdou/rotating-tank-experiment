clear all;
close all;
clc;

% this script generate a 2D grid in (X,Y) coordinate system for the
% simulation of psi-zeta equation (1.5-layer model for the beta-plume generated
% by a localized source):
% \partial/\partial t(\zeta-kd^2*\psi)=J(\zeta,\psi)-\beta\psi_x+Qf
% where Qf represents the mass source in terms of a Gaussian function of
% space.%The flow is confined in a cylindrical tank, where a solid barrier
% is placed inside to the west of a localzied source. The barrier is
% aligned in approximately the radial direction, with
WbOclock=1; %denoting the orientation (pointing to WbOclock o'clock)

aiv_im=imread(['for_grid_gen/',num2str(WbOclock),'oclock.tif']);

figure;
imshow(aiv_im);
hold on;
[px,py]=ginput(4); %left right top bot
xc=mean(px(1:2));yc=mean(py(3:4));
plot(xc,yc,'ys')
[sx,sy]=ginput(1); %source
plot(sx,sy,'ro','markersize',10)
[bx,by]=ginput(2); %western barrier
plot(bx,by,'m--','linewidth',5)
close;

nv2source= sx-xc - 1j*(sy-yc);
nv2source=nv2source./abs(nv2source);

nv2barrier= bx-xc - 1j*(by-yc);
dbar=110*abs(nv2barrier(1))/ abs(diff(px(1:2))); %cm distance from the tip of the barrier to the center

nv2barrier=nv2barrier./abs(nv2barrier);



% diameter of the cylindrical tank 
tank_dia = 110;
swthick = 5;%side wall thickness
Lx = tank_dia+swthick;% cm, the additional 5 cm is used to implement the boundary condtions
N = 500+1;% grid resolution in X and Y directoin
grid_filename = ['N',num2str(N-1),' BarrierOiren ',num2str(WbOclock,'%2.1f'),'.mat'];

i=2:N-1; %x-indices, increasing along x-aixs
j=i; %y-indices, increasing along y-axis
x=linspace(0,Lx,N);  %x-grid used for drawing pictures
[Y,X]=meshgrid(x,x); %horizonal grid
r=sqrt((X-Lx*0.5).^2+(Y-Lx*0.5).^2); %radial distance to the axis of rotation, it is used to identify the boundary of the domain
dx=Lx/(N-1);   %spatial resolution, cm
rdx=1./dx;


% find the solid points
Wall=r>=55;
figure;
imagesc(x,x,Wall');
ylabel('y');
xlabel('x');
set(gca,'ydir','normal')
hold on;
plot(Lx*0.5+[real(dbar*nv2barrier(1)) real(55*nv2barrier(2)) ],......
     Lx*0.5+[imag(dbar*nv2barrier(1)) imag(55*nv2barrier(2)) ],......
    'y-','linewidth',5);
plot(Lx*0.5,Lx*0.5,'ys')
% place the wall(a bar represneted an rectangle)
for n=1:4
    [wall_xy(n,1),wall_xy(n,2)]=ginput(1);
    text(wall_xy(n,1),wall_xy(n,2),num2str(n),'color','r')
end
patch(wall_xy(:,1), wall_xy(:,2), [0.5 0.5 0.5],'FaceAlpha', .1)
Wall=Wall|inpolygon(X,Y,wall_xy(:,1),wall_xy(:,2)) ;


%%%%%%%%%%define the boundary of the water domian such that psi(boundary>0)==0 maintains
boundary=zeros(N);
boundary(i,j) = ~Wall(i,j)&Wall(i+1,j)|~Wall(i,j)&Wall(i-1,j)|~Wall(i,j)&Wall(i,j+1)|~Wall(i,j)&Wall(i,j-1);
bp=find(boundary); %location of the boundary points in (X,Y) coordinate sys.
nbp=numel(bp);     %no. of boundary points
% inner and outside points used for boundary condition implementation.
bp_in_1=zeros(nbp,1); % -1 point inside the boundary(liquid points)
bp_in_3=zeros(nbp,1); % -3 point inside the boundary(liquid points)
bp_out_1=zeros(nbp,1);% +1 point outside the boundary(solid points/ ghost points)

% categorize the boundary points into 4 types:
% 1,side wall of the cylindrical tank: boundary = 1;
% 2,frontside of the western barrier facing the source: boundary = 2;
% 3,backside of the western barrier: boundary = 3;
% 4,tip side of the western barrier: boundary = 4;

figure;
imagesc(x,x,boundary');
ylabel('y');
xlabel('x');
set(gca,'ydir','normal',......
    'xlim',[min(wall_xy(:,1))-5,max(wall_xy(:,1))+5],......
    'ylim',[min(wall_xy(:,2))-5,max(wall_xy(:,2))+5],....
    'DataAspectRatio',[1 1 1]);
% choose mannually side2 points, which face the westward-propagating RWs
% radiated from the source region
for n=1:10
    [side2(n,1),side2(n,2)]=ginput(1);
    text(side2(n,1),side2(n,2),num2str(n),'color','y')
end
patch(side2(:,1), side2(:,2), [0.1 0.5 0.5],'FaceAlpha', .1)
boundary(bp(inpolygon(X(bp),Y(bp),side2(:,1),side2(:,2))))=2;

% choose mannually side3 points, which are on the backside the barrier
clf
imagesc(x,x,boundary');
ylabel('y');
xlabel('x');
set(gca,'ydir','normal',......
    'xlim',[min(wall_xy(:,1))-5,max(wall_xy(:,1))+5],......
    'ylim',[min(wall_xy(:,2))-5,max(wall_xy(:,2))+5],....
    'DataAspectRatio',[1 1 1]);
for n=1:10
    [side3(n,1),side3(n,2)]=ginput(1);
    text(side3(n,1),side3(n,2),num2str(n),'color','m')
end
patch(side3(:,1), side3(:,2), [0.5 0.5 0.1],'FaceAlpha', .1)
boundary(bp(inpolygon(X(bp),Y(bp),side3(:,1),side3(:,2))))=3;

% the left boundary points must located on the tip side of the barrier
boundary(r<45&boundary==1)=4;

% draw the boundaries
clf
imagesc(x,x,boundary');
ylabel('y,cm');
xlabel('x,cm');
set(gca,'ydir','normal',......
    'DataAspectRatio',[1 1 1],'xlim',[0 Lx],'ylim',[0 Lx]);
for s=2:4
    text( mean(X(boundary==s)),mean(Y(boundary==s)),['side ',num2str(s)],'color','w')
end
if all(boundary(bp)>0)
else
    error('some boundary points are not asigned')
end
title('4 sides');
text(Lx/2,5,'side 1','color','w')
% print(gcf,'grid and 4 sides','-dpng')










% %%%%%%%%%%%%%%%%% find the boundary, inside, and ghost points %%%%%%%%%%%%
ij=([find(~Wall&boundary==0) ]);  % inside points excluding the boundary points
Ninsider=numel(ij);


% Now let's locate the inner and outer points in the (X,Y) grid.
hw=2; %half window size used for identifying the inside and ghost/outside points
hl=5;% half length of the local tangential line
hs=3;%cm

figure('position',[0,0,600,600],'name','bp points');
axes('position',[0,0,1,1])
plot(X(ij),Y(ij),'marker','.','linestyle','none','color',[0.05 0.5 0.5 0.1]);
hold on;
plot(X(bp),Y(bp),'marker','s','linestyle','none','color','k');


for ii=1:nbp
    
    dist=(X(bp)-X(bp(ii))).^2 + (Y(bp)-Y(bp(ii))).^2;
    [~,Ind]=sort(dist);
    [ix,jy]=ind2sub([N,N],bp(Ind(1)));
    
    % a local tangential line obtained by linear regression
    xtang=X(bp(Ind(1:hl))) ;
    ytang=Y(bp(Ind(1:hl))) ;
    mx=mean(xtang);
    my=mean(ytang);        
%     plot(X(ix,jy), Y(ix,jy), 'gs','markerfacecolor','g');
    
    
    %   get the normal vector direction to the tangential line
    if boundary(bp(ii))==1 %circular side wall of the tank
        % normal vector in radial direction (pointing to the tank center)
        p(1)=( Lx*0.5-ytang(1) )/ (  Lx*0.5-xtang(1) );
        p(2)=ytang(1)-p(1)*xtang(1);
    else
        if all(xtang==xtang(1)) %tangential line is vertical
            p(2)=ytang(1);
            p(1)=0;        %slope for the normal vector to the tangantial line
        else    %tangential line is along sides :2,3,4
            [p,S] = polyfit(xtang,ytang,1);
            p(2)=ytang(1)+1./p(1)*xtang(1); %y intercept
            p(1)=-1./p(1);
        end
    end
    
    xf=linspace( mx-3,mx+3,10);
    yf=polyval(p,xf);
    plot(xf,yf,'color',[0.5 0.5 0.5 0.5]);% draw those normal vectors
    
    
    % unit vectors normal to the 4 solid bounadries, and they point inward towards the
    % liquid region
    if boundary(bp(ii))==1
        cosc=1/sqrt(1+p(1)^2)*sign( Lx*0.5-X(bp(Ind(1))) );
        sinc=abs(p(1)*cosc)*sign( Lx*0.5-Y(bp(Ind(1))) );
    elseif  boundary(bp(ii))==2 % frontside of the western barrier,facing the west-propagating RWs from the source
        if WbOclock*(WbOclock-3)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        elseif (WbOclock-3)*(WbOclock-6)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        elseif (WbOclock-6)*(WbOclock-9)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        elseif (WbOclock-9)*(WbOclock-12)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        else
            error('Orientation of the barrier must between [0,12] O''clock')
        end
    elseif  boundary(bp(ii))==3 % backside of the western barrier
        if WbOclock*(WbOclock-3)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        elseif (WbOclock-3)*(WbOclock-6)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        elseif (WbOclock-6)*(WbOclock-9)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        elseif (WbOclock-9)*(WbOclock-12)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        else
            error('Orientation of the barrier must between [0,12] O''clock')
        end
    elseif  boundary(bp(ii))==4  % tip side of the western barrier
        if WbOclock*(WbOclock-3)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        elseif (WbOclock-3)*(WbOclock-6)<=0
            cosc = -1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        elseif (WbOclock-6)*(WbOclock-9)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = abs(p(1)*cosc);
        elseif (WbOclock-9)*(WbOclock-12)<=0
            cosc = 1/sqrt(1+p(1)^2);
            sinc = -abs(p(1)*cosc);
        else
            error('Orientation of the barrier must between [0,12] O''clock')
        end
    else
        error('only 4 boundaries are needed.')
    end
    
    %find the inner points (1 dx away from the boundaries)
    vect= X(ix-hw:ix+hw,jy-hw:jy+hw)-X(ix,jy) +1j*((Y(ix-hw:ix+hw,jy-hw:jy+hw)-Y(ix,jy)));
    vect_n=vect./abs(vect);
    dd= abs(vect)*rdx;
    dd(dd==0)=nan;
    ang= real((cosc-sinc*1j)*vect_n); % max <=1
    ang(ang<=0)=nan;
    vect=abs(abs(dd-1)+3j*abs(ang-1)) ;
    [~,Ind2]=min(vect(:));
    [ni,nj]=ind2sub([2*hw+1,2*hw+1],Ind2);
    hold on
    % plot(X2(Ind2),Y2(Ind2),     'rs','markerfacecolor','r');
    % text(X2(Ind2),Y2(Ind2),num2str(ii))
    plot([ X(ix+ni-(hw+1),jy+nj-(hw+1))],.....
        [ Y(ix+ni-(hw+1),jy+nj-(hw+1))],.....
        'bo','markerfacecolor','b');
    xxx= sub2ind([N,N],ix+ni-(hw+1),jy+nj-(hw+1));
    if ~ismember(xxx,ij)
        %mannually locate the ghost point
        plot(X(ix,jy),Y(ix,jy),'cp','markerfacecolor','c')
        set(gca,'xlim', [X(ix,jy)-hs X(ix,jy)+hs],'ylim', [Y(ix,jy)-hs Y(ix,jy)+hs])
        set(gcf,'name','in 1 ')
        [px,py]=ginput(1);
        dist=(X-px).^2+(Y-py).^2;
        [~,xxx]=min( dist(:));
        bp_in_1(ii)=xxx;
        plot(X(xxx),Y(xxx),'bp','markerfacecolor','b');
    else
        bp_in_1(ii)=xxx;
    end
    
    %find the inner points (3 dx away from the boundaries)
    vect=abs(abs(dd-3)+3j*abs(ang-1)) ;
    [~,Ind2]=min(vect(:));
    [ni,nj]=ind2sub([2*hw+1,2*hw+1],Ind2);
    hold on
    % plot(X2(Ind2),Y2(Ind2),     'rs','markerfacecolor','r');
    % text(X2(Ind2),Y2(Ind2),num2str(ii))
    plot([ X(ix+ni-(hw+1),jy+nj-(hw+1))],.....
         [ Y(ix+ni-(hw+1),jy+nj-(hw+1))],.....
        'mo','markerfacecolor','m'); %inward -3 point
    xxx= sub2ind([N,N],ix+ni-(hw+1),jy+nj-(hw+1));
    if ~ismember(xxx,ij)
        plot(X(ix,jy),Y(ix,jy),'cp','markerfacecolor','c')
        set(gca,'xlim', [X(ix,jy)-hs X(ix,jy)+hs],'ylim', [Y(ix,jy)-hs Y(ix,jy)+hs])
        set(gcf,'name','in 3 ')
        [px,py]=ginput(1);
        dist=(X-px).^2+(Y-py).^2;
        [~,xxx]=min( dist(:));
        bp_in_3(ii)=xxx;
        plot(X(xxx),Y(xxx),'mp','markerfacecolor','m');
    else
        bp_in_3(ii)=xxx;
    end
        
    %find the outer/ghost points (1 dx away from the boundaries)
    ang= real((cosc-sinc*1j)*vect_n); % max <=1
    ang(ang>=0)=nan;
    vect=abs(abs(dd-1)+3j*abs(ang+1)) ;
    [~,Ind2]=min(vect(:));
    [ni,nj]=ind2sub([2*hw+1,2*hw+1],Ind2);
    hold on
    % plot(X2(Ind2),Y2(Ind2),     'rs','markerfacecolor','r');
    % text(X2(Ind2),Y2(Ind2),num2str(ii))
    plot([ X(ix+ni-(hw+1),jy+nj-(hw+1))],.....
        [ Y(ix+ni-(hw+1),jy+nj-(hw+1))],.....
        'ys','markerfacecolor','y'); %outward
    xxx=sub2ind([N,N],ix+ni-(hw+1),jy+nj-(hw+1));
    %  The ghost points cannot be coincide with the boundary points,or drop inside the
    % liquid area;
    % For the 1st case, psi(bp)=0 is maintained during time integration
    
    if ismember(xxx,ij)
        %mannually locate the ghost point
        plot(X(ix,jy),Y(ix,jy),'cp','markerfacecolor','c')
        set(gca,'xlim', [X(ix,jy)-hs X(ix,jy)+hs],'ylim', [Y(ix,jy)-hs Y(ix,jy)+hs])
        set(gcf,'name','out 1 ')
        [px,py]=ginput(1);
        dist=(X-px).^2+(Y-py).^2;
        [~,xxx]=min( dist(:));
        bp_out_1(ii)=xxx;
        plot(X(xxx),Y(xxx),'yp','markerfacecolor','y');
    else
        bp_out_1(ii)=xxx;
    end
    
end
set(gca,'xlim', [0 Lx],'ylim',[0 Lx]);








% all the points needed 
ij=([find(~Wall&boundary==0) ; bp ;bp_out_1]);  % all points need to be updated, all the others are treated as zero
NN=numel(ij);

% a double-check, the following 4 if must be checked before running the
% model
if sum(ismember(bp, ij))~=numel(bp);
    error('some boundary points are not included in ij')
end
if sum(ismember(bp_in_1, ij))~=numel(bp);
    error('some inside-1 points are not included in ij')
end
if sum(ismember(bp_in_3, ij))~=numel(bp);
    error('some inside-3 points are not included in ij')
end
if sum(ismember(bp_out_1, ij))~=numel(bp);
    error('some outside+1 points are not included in ij')
end


% %-----location of the neighbouring points for the inside points (ix, iy)
[ix,iy]=ind2sub( [N,N], ij(1:Ninsider));%care about the inside points

% right 
[~,Locb_r]  = ismember(sub2ind( [N,N],ix+1,iy), ij);
Locb_r(Locb_r==0) = NN+1; %zero is asigned and maintained there(where ij exclueds)

% left
[~,Locb_l]  = ismember(sub2ind( [N,N],ix-1,iy), ij);
Locb_l(Locb_l==0)=NN+1;

% upper
[~,Locb_u]  = ismember(sub2ind( [N,N],ix,iy+1), ij);
Locb_u(Locb_u==0)=NN+1;

% bottom
[~,Locb_b]  = ismember(sub2ind( [N,N],ix,iy-1), ij);
Locb_b(Locb_b==0)=NN+1;

% upper-right
[~,Locb_ur] = ismember(sub2ind( [N,N],ix+1,iy+1), ij);
Locb_ur(Locb_ur==0)=NN+1;

% upper-left
[~,Locb_ul] = ismember(sub2ind( [N,N],ix-1,iy+1), ij);
Locb_ul(Locb_ul==0)=NN+1;

% bottom-left
[~,Locb_bl] = ismember(sub2ind( [N,N],ix-1,iy-1), ij);
Locb_bl(Locb_bl==0)=NN+1;

% bottom-right
[~,Locb_br] = ismember(sub2ind( [N,N],ix+1,iy-1), ij);
Locb_br(Locb_br==0)=NN+1;

%  boundary points
[~,Loc_bp] = ismember(bp, ij);
% inner points dropping inside of the boundary(liquid region), with distance 1dx to
% the boundary
[~,Loc_bp_in1] = ismember(bp_in_1, ij);
% inner points dropping inside of the boundary(liquid region), with distance 3*dx to
% the boundary
[~,Loc_bp_in3] = ismember(bp_in_3, ij); 
% outter/ghost points dropping outside of the boundary (solid region), with distance 1*dx away
% from the boundary
[~,Loc_bp_out]= ismember(bp_out_1, ij); %bp_out_1 is not contained in ij


% if N^2>2^32
%     ij=int64(ij);
%     
%     Locb_b=int64(Locb_b);
%     Locb_u=int64(Locb_u);
%     Locb_l=int64(Locb_l);
%     Locb_r=int64(Locb_r);
%     
%     Locb_bl=int64(Locb_bl);
%     Locb_ul=int64(Locb_ul);
%     Locb_br=int64(Locb_br);
%     Locb_ur=int64(Locb_ur);
%     
% Loc_bp=int64(Loc_bp);
% Loc_bp_in1=int64(Loc_bp_in1);
% Loc_bp_in3=int64(Loc_bp_in3);
% Loc_bp_out=int64(Loc_bp_out);
% elseif N^2>2^16
%      ij=int32(ij);
%     
%     Locb_b=int32(Locb_b);
%     Locb_u=int32(Locb_u);
%     Locb_l=int32(Locb_l);
%     Locb_r=int32(Locb_r);
%     
%     Locb_bl=int32(Locb_bl);
%     Locb_ul=int32(Locb_ul);
%     Locb_br=int32(Locb_br);
%     Locb_ur=int32(Locb_ur);
%     
% Loc_bp=int32(Loc_bp);
% Loc_bp_in1=int32(Loc_bp_in1);
% Loc_bp_in3=int32(Loc_bp_in3);
% Loc_bp_out=int32(Loc_bp_out);
% 
% else
%       ij=int16(ij);
%     
%     Locb_b=int16(Locb_b);
%     Locb_u=int16(Locb_u);
%     Locb_l=int16(Locb_l);
%     Locb_r=int16(Locb_r);
%     
%     Locb_bl=int16(Locb_bl);
%     Locb_ul=int16(Locb_ul);
%     Locb_br=int16(Locb_br);
%     Locb_ur=int16(Locb_ur);
%     
% Loc_bp=int16(Loc_bp);
% Loc_bp_in1=int16(Loc_bp_in1);
% Loc_bp_in3=int16(Loc_bp_in3);
% Loc_bp_out=int16(Loc_bp_out);
% end
    


save(grid_filename,'Loc*','Lx','N','Ninsider','ij','WbOclock','Wall','nv2source','nv2barrier','dbar');










