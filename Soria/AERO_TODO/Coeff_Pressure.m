% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function [cp] = Coeff_Pressure (datos)
    
        dx = datos.L/datos.Vx;
        dy = datos.H/datos.Vy;
        Nx = datos.Nx;
        Ny = datos.Ny; 
        
    cp.ae = zeros(datos.Nx);
    cp.aw = zeros(datos.Nx);
    cp.an = zeros(datos.Nx);
    cp.as = zeros(datos.Nx);
    cp.ap = zeros(datos.Nx);
    cp.bp = zeros(datos.Nx);
    
    for i=2:Nx-1
        for j=2:Ny-1
            cp.ae(i,j) = dy/dx;
            cp.aw(i,j) = dy/dx;
            cp.an(i,j) = dx/dy;
            cp.as(i,j) = dx/dy; 
            cp.bp(i,j)= 0; 
        end
    end
    
    cp.ap = cp.ae + cp.aw + cp.an + cp.as;