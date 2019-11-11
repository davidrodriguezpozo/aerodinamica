function p = Analytic_Pressure (datos, C, time)
% This function ... wh
    for i = 1:datos.Nx
        for j=1:datos.Ny
        
        mu = datos.mu;
        rho = datos.rho;
        x = C.X(i);    
        y = C.Y(j);
        p(i,j) = -(exp(-8*pi^2*mu*time))^2 * rho * (cos(4*pi*x)+ cos(4*pi*y))/4;
           
        end 
    end