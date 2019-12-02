function [dt] = CalculateTimeStep (d, u, v)

    dx = d.L/d.Nx;
    
    Conv_dt = 0; Diff_dt = 0;
    
    umax = 0; vmax = 0;
    
    for i = 1:d.Nrow_sx
        for j = 1:d.Ncol_sx
            if u(i,j) > umax
                umax = u(i,j);
            end
        end
    end
    
    for i = 1:d.Nrow_sy
        for j = 1:d.Ncol_sy
            if v(i,j) > vmax
                vmax = v(i,j);
            end
        end
    end
    
    modVel = sqrt (umax*umax + vmax*vmax);
    
    Conv_dt = 0.35*dx/modVel;
    Diff_dt = 0.2*d.rho*dx*dx/d.mu; 
    
    dt = min(Conv_dt ,Diff_dt)*0.8; 
    