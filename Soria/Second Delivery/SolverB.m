%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION SolverB        %%%%%%%%%%%%%
%%% This function allows enables the programmer to      %%%
%%%   obtain the pressure field solution for a given    %%%
%%% velocity field for different sizes of the mesh Vx.  %%%
%%% As well as the maximum error for this mesh of       %%%
%%%                   discretization.                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function [u, v, p, step_time] = SolverB(datos, C, u_p, v_p)
    delta_t = 0;
    [nodal_mesh, num] = nodalmesh(datos.Vx,datos.Vy);
    [p, pseudo_p, step_time] = pressure (datos, u_p, v_p, nodal_mesh, delta_t);
    [p_gradX, p_gradY] = gradient_p(datos, pseudo_p, nodal_mesh); 
    p_gradX = haloupdate(p_gradX);
    p_gradY = haloupdate(p_gradY);
    u_p = haloupdate(u_p);
    v_p = haloupdate(v_p);

    %Corrected velocities
    u = u_p - p_gradX;
    v = v_p - p_gradY;

    divergence_field = divergencia_u(datos, u, v, nodal_mesh, num);
    % 
    sum(divergence_field)
    
    