function p = Analytic_Pressure (datos, C)
% This function ... wh
    for i = 1:datos.Nx
        for j=1:datos.Ny
        
        x = C.X(i);    
        y = C.Y(j);
        p(i,j) = - datos.F^2 * datos.rho * (cos(4*pi*x)+ cos(4*pi*y))/4;
           
        end 
    end