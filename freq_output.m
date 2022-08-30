


% if mod(Niter,Nnt)==0 %draw pictures
    clear slice*
    out_fname=['t ',num2str(T/(2*pi/OM),'%3.1f'),'r'];
    
    
    
    % figure('name','eta,cm','position',[10,10,1200,500])
%    h1=axes('Position',[0.0,0.005, 0.44 ,1100*0.45/500]);
%    slice_p=zeros(N);
%    slice_p(ij(1:Ninsider))=psi(1:Ninsider)*f_Coriolis/g;
%    v=slice_p(3:end,2:end-1)-slice_p(1:end-2,2:end-1);
%    v=v*g/f_Coriolis;
%    u=slice_p(2:end-1,3:end)-slice_p(2:end-1,1:end-2);
%    u=u*g/f_Coriolis;
%    u(u.*v==0)=nan;
    
%    imagesc(x,x,slice_p','AlphaData',~Wall')   
    %         ylabel('y, cm');
    %         xlabel('x, cm');
%    hold on;
%    quiver(x(2:nv:end-1),x(2:nv:end-1),-u(2:nv:end-1,2:nv:end-1)'*uv_scale,v(2:nv:end-1,2:nv:end-1)'*uv_scale,0,'color','k')
%    quiver(96,10,1*uv_scale,0,0,'color','w','linewidth',1,'maxheadsize',18)
%    text(90,4,'1 cm/s','fontsize',18,'fontweight','bold','color','w')    
%     text(5,108,['eta, cm. ',num2str(round(T)),'s'],'color','w','fontsize',18);
%    set(gca,'ydir','normal','color','k','xtick',[],'ytick',[]);
%    hb=colorbar;hb.Location='eastoutside';hb.Color='k';
%    hb.Position=[0.445,0.05,0.01,1-0.05*2];
%    hb.FontSize=12;
%  text(1,111,['\gamma=',num2str(max(Ek_gamma),'%6.4f'),'s^{-1}'],'fontsize',20,'color','w');
%  text(1,101,['\nu=',num2str(Ah,'%6.4f'),'cm^2/s'],'fontsize',20,'color','w');
%  text(1,5,['Q=',num2str(Amp,'%6.4f'),'s^{-2}'],'fontsize',20,'color','w');

    
    
 %   h2=axes('Position',[0.507,0.005,0.44 ,1100*0.45/500]);
    h2b=axes('Position',[0,0,900/1000 ,1]);
    slice_z=zeros(N);
    slice_z(ij(1:Ninsider))=zeta(1:Ninsider)/f_Coriolis ;
    %         figure
    imagesc(x,x,slice_z' ,'AlphaData',(~Wall )');
    %     hb=colorbar;hb.Location='eastoutside';hb.Color='k';
    
%         title(['Ro=\zeta/f_0 ',subtitle]);
    caxis([-1 1]*0.2);
    %         xlabel('x,cm');
    %         ylabel('y,cm');
    set(gca,'ydir','normal','color','k','xtick',[],'ytick',[]);
%         text(5,108,['Ro '],'color','w','fontsize',18);
    hb=colorbar;
%    hb.Location='eastoutside';
    hb.Color='w';
   % hb.Position=[0.951,0.05,0.01,1-0.05*2];
    hb.Position=[0.89,0.05,0.01,0.9];
    hb.FontSize=12;
  text(1,111,[num2str(T/(2*pi/OM),'%3.0f'),' rot.'],'fontsize',20,'color','w','fontname','serif');
  text(99,109,['\zeta/f_0'],'fontsize',20,'color','w','fontname','serif');
    
    
    
    
    print(gcf,[dir2cf,case_folder,'/',out_fname,'.png'],'-dpng');        
%    save([dir2cf,case_folder,'/',out_fname],'psi');
    
      frame = getframe(gcf);
 writeVideo(v,frame);
    clf;

