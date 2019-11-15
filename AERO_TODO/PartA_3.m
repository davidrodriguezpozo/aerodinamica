% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function PartA_3

%set(groot, 'DefaultTextInterpreter','latex');
%set(groot, 'Defaultaxesticklabelinterpreter','latex');
%set(groot, 'DefaultLegendInterpreter','latex');

    clear all 
    clc
    close all 

    Vx = [3 5 10 20 30 50 100]; 
    Vy = [3 5 10 20 30 50 100]; 
    Re = 1000;
    
    for p=1:length(Vx)
    
    datos = INPUT(Vx(p),Vy(p),Re);
    datos.F = 1;
    C = meshes (datos, datos.malla);  
    
    u = zeros(datos.Nx, datos.Ny);
    v = zeros(datos.Nx, datos.Ny);
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1

            x1 = C.stagX_x(i,j);
            y1 = C.stagX_y(i,j);

            u(i,j) = datos.F*cos(2*pi*x1)*sin(2*pi*y1);

            x2 = C.stagY_x(i,j);
            y2 = C.stagY_y(i,j);

            v(i,j) = -datos.F*cos(2*pi*y2)*sin(2*pi*x2);

        end
    end

    u = haloupdate(u);
    v = haloupdate(v);
    
    [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);
    
    [conv_u_analytic diff_u_analytic] = Analytic (datos, C, 'x');
    [conv_v_analytic diff_v_analytic] = Analytic (datos, C, 'y');
    
    error_du(p) = ERROR(diff_u, diff_u_analytic);
    error_dv(p) = ERROR(diff_v, diff_v_analytic);
    error_cu(p) = ERROR(conv_u, conv_u_analytic);
    error_cv(p) = ERROR(conv_v, conv_v_analytic);
    
    end
    
    PlotsA (Vx,error_du,error_dv,error_cu,error_cv);