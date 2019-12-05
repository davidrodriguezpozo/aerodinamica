%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION SolverC        %%%%%%%%%%%%%
%%%   This function allows the programmer to obt        %%%
%%%   the pressure field solution and velocity          %%%
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

function [u, v, P, P_an, P_num, u_an, u_num, v_an, v_num, acu_time] =  SolverC(datos, C, ii, jj)

[u, v, R_u, R_v, u_p, v_p, u_prev, v_prev] = CondicionesIniciales(datos, C);
final_time = 10;
P = zeros(datos.Nx);
P_prev = zeros(datos.Nx);
dif_P = 1;
dif_U = 1;
dif_V = 1;

Nx = datos.Nx;
Ny = datos.Ny;

delta_t = 0.01;
delta_t = TimeStep (datos,u,v,delta_t);
time = 0;
R_uant = zeros(Nx,Nx);
R_vant = zeros(Nx,Nx);
index = 1;

while time <= final_time
    
    [nodal_mesh, num] = nodalmesh(datos.Vx,datos.Vy);
    [conv_u, diff_u, conv_v, diff_v] = Numerical (datos, C, u, v);
    [u_p, v_p, R_uant, R_vant] = predicted(u,v,delta_t,datos, conv_u, diff_u, conv_v, diff_v, R_uant, R_vant,R_u,R_v);
    
%     
%     R_u_prev = R_u;
%     R_v_prev = R_v;
    
    u_prev = u;
    v_prev = v;
    
    
    u_p = haloupdate(u_p);
    v_p = haloupdate(v_p);
    
    dx = datos.L/datos.Vx;
    dy = datos.H/datos.Vy;
    rho = datos.rho;
    mu = datos.mu;
    
%     [P, pseudo_p, step_time] = pressure (datos, u_p, v_p, nodal_mesh, delta_t);
%     [p_gradX, p_gradY] = gradient_p(datos, pseudo_p, nodal_mesh); 
    P = pressureC(datos, u, v, u_p, v_p, delta_t,P);
    [p_gradX, p_gradY] = gradient_pC(datos, P, nodal_mesh,delta_t);
    p_gradX = haloupdate(p_gradX);
    p_gradY = haloupdate(p_gradY);
%     u_p = haloupdate(u_p);
%     v_p = haloupdate(v_p);

    %Corrected velocities
    u = u_p - p_gradX;
    v = v_p - p_gradY;
    
    p_analytic = Analytic_Pressure (datos, C, time);
    rho = datos.rho;
    mu = datos.mu;
    
    u = haloupdate(u);
    v = haloupdate(v);

    u_num(index)=u(ii,jj);
    v_num(index)=v(ii,jj);
    P_an(index)= -(exp(-8*pi^2*mu*time))^2*(rho*(cos(4*pi*C.Coll_x(ii,jj))+cos(4*pi*C.Coll_y(ii,jj))))/4;
    P_num(index)=P(ii,jj);
    u_an(index) = datos.F*cos(2*pi*C.stagX_x(ii,jj))*sin(2*pi*C.stagX_y(ii,jj));
    v_an(index) = -datos.F*cos(2*pi*C.stagY_y(ii,jj))*sin(2*pi*C.stagY_x(ii,jj));
    acu_time(index) = time;
    index = index+1;

    dif_U = 0;
    dif_V = 0;
    %delta_t = 0.01;
    time = time + delta_t;
    
    datos.F = exp(-8*pi^2*datos.mu*time);
    delta_t = TimeStep (datos,u,v,delta_t)

end