%% DRIVEN CAVITY PROBLEM 

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

% Execute function main 

function PartD_3

clear all
clc
close all

d = INPUT;

C = MESH (d);

S = CondicionesIniciales(d,C);

[cp] = Coefficients_Pressure(d);

S = Solve_NS(d,C, S, cp);

Plots(C,S.U,S.V,S.P);
