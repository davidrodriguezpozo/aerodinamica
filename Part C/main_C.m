function main_C

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 

%                                 El Codi hauria de fer aixo 
% 
%  1. Timestep for stability
%  2. u_n_1 i u_n conegudes --> R_n_1 i R_n
%  3. Evaluar u_p = u_n + delta_t(3/2 * R_n - 1/2 * R_n_1
%         R = - conv(u)/V + diff(u)/V
%  4. Resoldre l'equaci? de Poisson A * pseudo_p = Div(u_p)
%  5. Obtenir pseudo_p
%  6. Obtenir u_n1 --> u_n1 = u_p - grad(pseudo_p) 
%  7. timestep + 1 


Vx = 10;

datos = Input(Vx,Vx);
C = meshes (datos, datos.malla);

time = datos.dt_inicial;
datos.F = exp( -8*pi^2 * datos.mu * time);
[u v] = Velocitats (datos, C);

while time < datos.time_final || difP < datos.delta

% Fem Update del temps¡
datos.F = exp( -8*pi^2 * datos.mu * time);
S.R_u_prev = S.R_u;
S.R_v_prev = S.R_v;

u_prev = u;
v_prev = v;

% PERFIL VELOCITATS
[u v] = Velocitats (datos, C);

[S.R_u S.R_v S.R_u_prev S.R_v_prev] = ComputeR (datos, C, u, v);

time = 1000;

end



