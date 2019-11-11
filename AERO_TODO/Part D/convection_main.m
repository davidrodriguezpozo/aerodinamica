function cu = convection_main (u,v,L,H)
% The inputs: the velocity field (u,v), and the initial values (L,H).
% The outputs: the convection.
        Nx = size(u,1);
        Ny = size(u,2);
        
        % The Alpha X and Alpha Y are computed with the initial values initialized in the main. 
        dx = L/Nx; dy = H/Ny;
        
% The double loop (for each direction) starts at i=2 and ends at i=N-1, which is the next-to-last node. 
    for i = 2:Nx-1
        for j=2:Ny-1

        % Using CDS (Central Difference Scheme) the components of the velocity in each surface of the VC are computed. 
        ue = (u(i+1,j)+u(i,j))/2;
        uw = (u(i+1,j)+u(i,j))/2;
        vn = (v(i,j+1)+u(i,j))/2;
        vs = (v(i,j-1)+u(i,j))/2;
        un = (u(i,j+1)+u(i,j))/2;
        us = (u(i,j-1)+u(i,j))/2;
   
        % Mass Fluxes using Alpha Y and Alpha X (the distance incrementals) 
        Fe = ue * dy;
        Fw = uw * dy;
        Fn = vn * dx;
        Fs = vs * dx;
    
        % Finally, the Convection term (for i-1 positions)
   
        cu(i,j) = Fe*ue - Fw*uw + Fn*un - Fs*us;
        end
    end