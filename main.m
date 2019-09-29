function main

clear all 
clc
close all 

global C
global datos
global S

datos = INPUT;
C = MESH (datos);
S = CondicionesIniciales(datos,C);
S = Solve_NS(datos,C);


function datos = INPUT 

    datos.Nx = 5;
    datos.Ny = 5;
        datos.Ncol = datos.Nx + 2; % Numero de nodos en x-direction
        datos.Nrow = datos.Ny + 2; % Numero de nodos en x-direction
        datos.Ncol_sx = datos.Ncol - 1;
        datos.Nrow_sx = datos.Nrow - 2;
        datos.Ncol_sy = datos.Ncol - 2;
        datos.Nrow_sy = datos.Nrow - 1;
    datos.L = 1;
    datos.H = 1;
    datos.uniform = true;
    datos.gamma = 10;
    datos.tfinal = 100;
    datos.dt_inicial = 0.01;
    