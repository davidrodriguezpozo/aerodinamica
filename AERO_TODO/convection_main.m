function [cu cv] = convection_main (u,v,L,H,C)
    

        Nx = size(u,1);
        Ny = size(u,2);
        cu = zeros(Nx,Ny);
        cv = zeros(Nx,Ny);
        
        Vx = Nx -2;
        Vy = Ny -2;
        
        dx = L/Vx; dy = H/Vy;

    for i = 2:Nx-1
        for j=2:Ny-1

        % CDS : Central Difference Scheme
        ue = (u(i,j)+u(i,j+1))/2;
        uw = (u(i,j)+u(i,j-1))/2;
        un = (u(i,j)+u(i+1,j))/2;
        us = (u(i,j)+u(i-1,j))/2;
        
        ve = (v(i,j)+v(i,j+1))/2;
        vw = (v(i,j)+v(i,j-1))/2;
        vn = (v(i,j)+v(i+1,j))/2;
        vs = (v(i,j)+v(i-1,j))/2;
   
        % Mass Fluxes 
        
        rho = 1;
        
        dx1 = dx/2;
        dy1 = dx/2;

        Fe_u = 0.5*rho*( u(i,j+1) + u(i,j))*(dy);
        Fw_u = 0.5*rho*(u(i,j-1)+u(i,j))*(dy);
        Fn_u = rho*v(i,j)*(dx1)+rho*v(i,j+1)*(dx1);
        Fs_u = rho*v(i-1,j)*(dx1)+rho*v(i-1,j+1)*(dx1);
        
        Fe_v =rho*u(i,j)*(dy1)+rho*u(i+1,j)*(dy1);
        Fw_v = rho*u(i,j-1)*(dy1)+rho*u(i+1,j-1)*(dy1);
        Fn_v = 0.5*(rho*v(i+1,j)+rho*v(i,j))*(dx);
        Fs_v = 0.5*(rho*v(i-1,j)+rho*v(i,j))*(dx);
        
        Fe = ue * dy;
        Fw = uw * dy;
        Fn = vn * dx;
        Fs = vs * dx;
    
        % Convection term
   
        cu(i,j) = Fe_u*ue - Fw_u*uw + Fn_u*un - Fs_u*us;
        cv(i,j) = Fe_v*ve - Fw_v*vw + Fn_v*vn - Fs_v*vs;
        
        end
    end