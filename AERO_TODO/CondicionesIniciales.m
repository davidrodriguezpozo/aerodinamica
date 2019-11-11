%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      FUNCTION CONDICIONES INICIALES     %%%%%%%
%%% This function sumarizes all the initial conditions  %%%
%%% in order have them ordered. It uses the structure   %%%
%%% "datos", which has the inputs, and then inisializes %%%
%%% the u and v velocitiy in a pre-determined velocity  %%%
%%% field. The haloupdated function is executed.        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod


    function [u v R_u R_v u_p v_p u_prev v_prev] = CondicionesIniciales(datos, C)
        
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

    R_u = zeros(datos.Nx);
    R_v = zeros(datos.Nx);

    conv_u = zeros(datos.Nx); 
    diff_u = zeros(datos.Nx);
    conv_v = zeros(datos.Nx);
    diff_v = zeros(datos.Nx);

    u_p = zeros(datos.Nx);
    v_p = zeros(datos.Nx);
    u_prev = zeros(datos.Nx);
    v_prev = zeros(datos.Nx);
