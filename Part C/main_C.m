function main_C

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 

Vx = 10;

datos = Input(Vx,Vx);
C = meshes (datos, datos.malla);

time = datos.dt_inicial;
datos.F = exp( -8*pi^2 * datos.mu * time);
[u v] = Velocitats (datos, C);

while time < datos.time_final || difP < datos.delta

% Fem Update del temps
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



