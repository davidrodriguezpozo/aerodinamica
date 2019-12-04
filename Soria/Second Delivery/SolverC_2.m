%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION SolverC        %%%%%%%%%%%%%
%%%   This function allows enables the programmer to    %%%
%%%   obtain the pressure field solution and velocity   %%%
%%%   field for each instant of time until a stabilty   %%%
%%%   is achieved or when final time is reached.        %%%
%%%   The program also computes the numerical and       %%% 
%%%   analytical value of a point in the domain.        %%%
%%%   Notice that this program is running along time.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function [u v P P_an P_num u_an u_num v_an v_num acu_time] =  SolverC(datos, C, ii, jj)

[u v R_u R_v u_p v_p u_prev v_prev] = CondicionesIniciales(datos, C);
final_time = 3;
[cp] = Coeff_Pressure (datos);

P = zeros(datos.Nx);
P_prev = zeros(datos.Nx);
dif_P = 1;
dif_U = 1;
dif_V = 1;
delta_V = 1e-7;

Nx = datos.Nx;
Ny = datos.Ny;

Td = 0.5*(datos.L/datos.Vx)*(datos.H/datos.Vy) / datos.mu;
delta_t = TimeStep(datos,u,v);
time = delta_T;

    
index = 1;

while(time<2)
        
    %Numerical Solution for the convective and diffusive terms:
    [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);

    %Predicted velocity calculation:
    [u_p, v_p, R_uant, R_vant] = predicted(u,v,delta_t,datos, conv_u, diff_u, conv_v, diff_v, R_uant, R_vant,R_u,R_v);
    u_p = haloupdate(u_p);
    v_p = haloupdate(v_p);
    
    
    [P,P0]=pressure_solver(P, P0, rho, up, vp, x_point, y_point, ux, vy, delta_t, volumes_x, volumes_y);
    [u, v]=velocity_calculation(P, delta_t, rho, vp, up, x_point, y_point, ux, vy, volumes_x, volumes_y, u, v);
    u=halo_update(u);
    v=halo_update(v);
    [t_vector, u_vector, v_vector, U_anal, V_anal, P_anal, P_vec]=comparison_analytical_solution(time, delta_t, index,...
    t_vector, u_vector, v_vector, U_anal, V_anal, P_anal, P_vec, P, u, v, ux, vy, mu, rho, x_point, y_point);
    [delta_t]=delta_t_calculation(x_point, y_point, u, v, volumes_x, volumes_y, delta_td, delta_t);     
    [max_error]=error_velocity_calculation(ux, vy, x_point, y_point, volumes_x, volumes_y, u, v, time, rho, mu);
    
    if(max_error>max_error_timestep) 
        max_error_timestep=max_error;
    end
    
    time=time+delta_t;
    index=index+1;
    
end

end