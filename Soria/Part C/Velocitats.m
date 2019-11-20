function [u v ] = Velocitats (datos, C)

for i = 1:datos.Nx
    for j = 1:datos.Ny
        
        x = C.stagX_x(i,j);
        y = C.stagX_y(i,j);
        
        u(i,j) = datos.F*cos(2*pi*x)*sin(2*pi*y);
        
        x = C.stagY_x(i,j);
        y = C.stagY_y(i,j);
        
        v(i,j) = -datos.F*cos(2*pi*y)*sin(2*pi*x);
        
    end
end