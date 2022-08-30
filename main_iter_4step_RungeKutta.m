
drawnow;
if Niter==0
   movie_name = [case_folder,' ',num2str(T_tot),'s']   ;
  % Hfigure=  figure('position',[100,10,1100,500],'InvertHardcopy','off','color','w');
   Hfigure=  figure('position',[100,10,500,500*900/1000],'InvertHardcopy','off','color','k');
  v=VideoWriter([dir2cf,case_folder,'/',movie_name]);
  v.FrameRate=30;
 open(v);


    if ~restart %the model starts from an initial field without motion
        %  such that the only term remains in the governning eq. is due to the external forcing
        % Euler forward for only the 1st time step.
        q1=dt*Qforcing;
        set(gcf,'Name',['New Run from T=',num2str(T,'%5.0f'),' seconds'])
    else
        set(gcf,'Name',['Restart Run from T=',num2str(T,'%5.0f'),' seconds'])
    %      if restart %starts from the last moment of an exit run
   clear q1
   load([dir2cf,case_folder,'/',fls(end).name])
   var_list= who('-file',[dir2cf,case_folder,'/',fls(end).name]);    

   if ~ismember('q1',var_list) %q1(used for restart) is not saved for some reason 
q1(1:Ninsider,1)=(psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs-kds*psi(1:Ninsider);
   end


    end
end






if no_slip
    
    while T<=T_tot %run the model for 30 seconds
  
     
        %  ------------------------ k1 for Runge-Kutta scheme--------------------
        
        %   Successive Over Relaxation method used to solve for psi in
        %      (\nabla^2 -kd^2) \psi = f*dt....
        %   ref:http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture11.pdf
        %    1:NN elements in psi are just for the inner points, excluding the
        %    boundary points, bp, and those dropping outside
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q1(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp_out )=psi( Loc_bp_in1 )*OrszagIsraeli74_22o13+ ....
        psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13;
        
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
        zeta( Loc_bp ) =  psi( Loc_bp_in1 )*OrszagIsraeli74_35o13 + ....% Orszag and Israeli 1974
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13  ;
        zeta( Loc_bp ) =    zeta( Loc_bp ) *rdxs*Ek_lateral;
        
  %  -------------------------beta term------------------------------------
        vr_beta = - coscita.*(psi(Locb_u)-psi(Locb_b)) + sincita.*(psi(Locb_r)-psi(Locb_l));
        vr_beta = vr_beta*rdx*0.5;
        vr_beta = beta_T.*vr_beta;    

        %         ------------------------ Ekman friction + Laplacian diffusion----
        Fric_Diff= -zeta(1:Ninsider).*Ek_gamma +.....Fric
            (zeta(Locb_r)+zeta(Locb_l)+zeta(Locb_u)+zeta(Locb_b)-4*zeta(1:Ninsider))*rdxs_Ah;     ......horizontal diffusion
        
        Jac=0;
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            Jac=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            Jac=Coe_Jac*Jac;
        end
        
 
        
        %  -----------------------k2 for Runge-Kutta scheme------------------------
        q2=q1+Jac*dtby2;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp_out )=psi( Loc_bp_in1 )*OrszagIsraeli74_22o13+ ....
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13;
        
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
%       zeta(1:Ninsider) = zeta(1:Ninsider)./(1+Ek_gamma*dtby2);
        zeta( Loc_bp ) =  psi( Loc_bp_in1 )*OrszagIsraeli74_35o13 + ....% Orszag and Israeli 1974
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13  ;
        zeta( Loc_bp ) =    zeta( Loc_bp ) *rdxs*Ek_lateral;
        
        q2=0;%jac2
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+2*q2;
        
        
        
        
        %  -----------------------k3 for Runge-Kutta scheme------------------------
        q2=q1+q2*dtby2;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp_out )=psi( Loc_bp_in1 )*OrszagIsraeli74_22o13+ ....
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13;
         
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
%       zeta(1:Ninsider) = zeta(1:Ninsider)./(1+Ek_gamma*dtby2);
        zeta( Loc_bp ) =  psi( Loc_bp_in1 )*OrszagIsraeli74_35o13 + ....% Orszag and Israeli 1974
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13  ;
        zeta( Loc_bp ) =    zeta( Loc_bp ) *rdxs*Ek_lateral;
        
        q2=0;%jac3
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+2*q2;
        
        
        %  -----------------------k4 for Runge-Kutta scheme------------------------
        q2=q1+q2*dt;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp_out )=psi( Loc_bp_in1 )*OrszagIsraeli74_22o13+ ....
            psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13;
        
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
%         zeta(1:Ninsider) = zeta(1:Ninsider)./(1+Ek_gamma*dt);
        zeta( Loc_bp ) =  psi( Loc_bp_in1 )*OrszagIsraeli74_35o13 + ....% Orszag and Israeli 1974
                          psi( Loc_bp_in3 )*OrszagIsraeli74_m1o13  ;
        zeta( Loc_bp ) =    zeta( Loc_bp ) *rdxs*Ek_lateral;
        
        q2=0;%jac4
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+q2;
        
        q2 = q1+dtby6*Jac+......nonlinear adv. Runge-Kutta scheme
            dt*(.......
            Qforcing +...... external forcing
            vr_beta +....    beta_T*vr
             Fric_Diff .....Laplacian diff, Ekman friction        
            ) ; 
        q1=q2;
        
        
        T=T+dt;
        Niter=Niter+1;
        
        if mod(Niter,Nnt)==0
            if any(isnan(psi(:)))
                disp(['model overflows at ',num2str(T,'%5.3f'),'s']);
                overflow=true;
                return
            else
                overflow=false;
            
            end
                    freq_output;
            disp([' T=',num2str(T,'%5.3f'),'s'])
%             set(gcf,'Name',['T=',num2str(T,'%5.0f'),' seconds'])
        end
    end
    
    
    if ~overflow
        out_fname=['t ',num2str(round(T)),'s'];
        save([dir2cf,'/',case_folder,'/',out_fname],'psi','q1');
        disp([dir2cf,'/',case_folder, '/',out_fname,' saved for restart!'])
    else
        disp(['this case overflows, no restart file will be created.'])
    end

    
    
else  %free slip b.c.
    
    
    
    
    
    
    
    while T<=T_tot %run the model for 30 seconds
        if T==0 %the model starts from an initial field without motion
            %  such that the only term remains in the governning eq. is due to the external forcing
            % Euler forward for only the 1st time step.
            q1=dt*Qforcing;
            
            Hfigure=  figure('position',[100,10,1100,500],'InvertHardcopy','off','Name','bt model','color','w');
            % axis tight manual
        end
        
        %  ------------------------ k1 for Runge-Kutta scheme--------------------
        
        %   Successive Over Relaxation method used to solve for psi in
        %      (\nabla^2 -kd^2) \psi = f*dt....
        %   ref:http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture11.pdf
        %    1:NN elements in psi are just for the inner points, excluding the
        %    boundary points, bp, and those dropping outside
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q1(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp )=0;
        
        zeta(1:Ninsider) = ( psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider) )*rdxs;         
        zeta( Loc_bp ) =0; %free slip condition
        
       %  -------------------------beta term------------------------------------
        vr_beta = - coscita.*(psi(Locb_u)-psi(Locb_b)) + sincita.*(psi(Locb_r)-psi(Locb_l));
        vr_beta = vr_beta*rdx*0.5;
        vr_beta = beta_T.*vr_beta;
   %         ------------------------ Ekman friction + Laplacian diffusion----
        Fric_Diff= -zeta(1:Ninsider).*Ek_gamma +.....Fric
            (zeta(Locb_r)+zeta(Locb_l)+zeta(Locb_u)+zeta(Locb_b)-4*zeta(1:Ninsider))*rdxs_Ah;     ......horizontal diffusion        
        
        
        Jac=0;
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            Jac=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            Jac=Coe_Jac*Jac;
        end
        

        
        %  -----------------------k2 for Runge-Kutta scheme------------------------
        q2=q1+Jac*dtby2;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp )=0;
        
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
        zeta( Loc_bp ) =0; %free slip condition
        
        q2=0;%jac2
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+2*q2;
        
        
        
        
        %  -----------------------k3 for Runge-Kutta scheme------------------------
        q2=q1+q2*dtby2;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp )=0;
        
        zeta(1:Ninsider) = ( psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider) )*rdxs;
        zeta( Loc_bp ) =0; %free slip condition
        
        q2=0;%jac3
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+2*q2;
        
        
        %  -----------------------k4 for Runge-Kutta scheme------------------------
        q2=q1+q2*dt;
        for n=1:50
            for ii=1:Ninsider % only for the inside points
                psi(ii)=  Coe_SOR*(.....
                    psi( Locb_l(ii)) + psi( Locb_r(ii)) + psi( Locb_u(ii)) + psi( Locb_b(ii)) .....
                    -dxs*q2(ii).....
                    )......
                    +(1-lambda_SOR)*psi(ii);
            end
        end
        psi( Loc_bp )=0;
        
        zeta(1:Ninsider) = (psi(Locb_r)+psi(Locb_l)+psi(Locb_u)+psi(Locb_b)-4*psi(1:Ninsider))*rdxs;
        zeta( Loc_bp ) = 0; %free slip condition
        
        q2=0;%jac4
        if ~linear
            % Arakawa's Jacobe(1966 eq(46)) is used for the nonlinear term J(\zeta,\psi)
            q2=  ( zeta(Locb_r)-zeta(Locb_l) ).*(psi(Locb_u)-psi(Locb_b))......
                - ( zeta(Locb_u)-zeta(Locb_b) ).*(psi(Locb_r)-psi(Locb_l))......
                +zeta(Locb_r).*(psi(Locb_ur)-psi(Locb_br))......
                -zeta(Locb_l).*(psi(Locb_ul)-psi(Locb_bl))......
                -zeta(Locb_u).*(psi(Locb_ur)-psi(Locb_ul))......
                +zeta(Locb_b).*(psi(Locb_br)-psi(Locb_bl))......
                +zeta(Locb_ur).*(psi(Locb_u)-psi(Locb_r))......
                -zeta(Locb_bl).*(psi(Locb_l)-psi(Locb_b))......
                -zeta(Locb_ul).*(psi(Locb_u)-psi(Locb_l))......
                +zeta(Locb_br).*(psi(Locb_r)-psi(Locb_b));
            q2=Coe_Jac*q2;
        end
        Jac=Jac+q2;
        
   
           q2 = q1+dtby6*Jac+......nonlinear adv. Runge-Kutta scheme
            dt*(.......
            Qforcing +...... external forcing
            vr_beta +....    beta_T*vr
             Fric_Diff .....Laplacian diff, Ekman friction        
            ) ; 
        q1=q2;
        
        T=T+dt;
        Niter=Niter+1;
        
        if mod(Niter,Nnt)==0
            if any(isnan(psi(:)))
                error(['model overflows at ',num2str(T,'%5.3f'),'s'])
            end            
            freq_output;
            disp([case_folder,' T=',num2str(round(T)),'s'])
        end
        
    end
    
end
close(v)
