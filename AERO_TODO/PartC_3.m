% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function PartC_3

%set(groot, 'DefaultTextInterpreter','latex');
%set(groot, 'Defaultaxesticklabelinterpreter','latex');
%set(groot, 'DefaultLegendInterpreter','latex');

clear all 
clc
close all 


Vx = 20; 
Vy = 20;
Re = 100;

datos = INPUT(Vx,Vy,Re);
C = meshes (datos, datos.malla);

datos.F = 1;

ii = 3; jj = 3;

[u v P P_an P_num u_an u_num v_an v_num acu_time] = SolverC(datos, C, ii, jj);

plotsC (u, v, P, P_an, P_num, u_an, u_num, v_an, v_num, acu_time, C, ii, jj, datos);


