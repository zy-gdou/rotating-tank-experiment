       
        % % the last element in psi is always 0, it is used for the boundary value(in case of the free slip b.c )
        % and is used for other useful points in case of the non-slip b.c.
        
%         initialzie matrices
        NN=numel(ij);
        psi=zeros(NN+1,1);
        zeta=zeros(NN+1,1);
        
        vr_beta=zeros(Ninsider,1); %beta_T.*vr
        q1=zeros(Ninsider,1);
        q2=zeros(Ninsider,1);
        Jac=zeros(Ninsider,1); %jac0
