%% NAVIER-STOKES COMPUTATIONA
% Code developed by:
% - Sergi Martínez Castellarnau
% - Carlos Pérez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

% Execute function main or main.m

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
    
    

    function C = MESH (datos)
    
    fprintf('-( 1 )- Calculando mallas...\n'); 
        
    % Preallocation of variables
    
    C.X = zeros (datos.Nrow,datos.Ncol);
    C.Y = zeros (datos.Nrow,datos.Ncol);    
    C.dx = zeros (datos.Nrow,1);
    C.dy = zeros (datos.Ncol,1);
    
    C.stagX_x = zeros (datos.Nrow_sx,datos.Ncol_sx);
    C.stagX_y = zeros (datos.Nrow_sx,datos.Ncol_sx);
    
    C.stagY_x = zeros (datos.Nrow_sy,datos.Ncol_sy);
    C.stagY_y = zeros (datos.Nrow_sy,datos.Ncol_sy);
    
    [ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( datos );

    [ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( datos, C );
    
    [ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( datos, C );
    
    figure; hold on;
    plot(C.Coll_x , C.Coll_y,'k.','MarkerSize', 25) ; 
    plot(C.stagX_x , C.stagX_y,'b.') ;
    plot(C.stagY_x , C.stagY_y,'g.') ;
    ylim ([-0.2 1.2])
    xlim ([-0.2 1.2])
    
    %legend('Collocated','','X-Staggered','','Y-staggered','');
    
    fprintf('---> Mallas calculada\n')      
    

function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( datos )
    
    % Delta X
    for i = 1:datos.Ncol 
            dx(i) = spacing_no_uniform( datos.gamma, datos.L, i, datos.Nx, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.L, i-1, datos.Nx, datos.uniform );
    end 
    
    % Delta Y
    for j = 1:datos.Nrow
            dy(j) = spacing_no_uniform( datos.gamma, datos.H, j, datos.Ny, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.H, j-1, datos.Ny, datos.uniform );
    end 
    
    % COLLOCATED MESHES
    
    X (1) = 0; X(2) = dx(1)/2;
    for i = 3:datos.Ncol-1
        X(i) = X(i-1) + dx(i-1);
    end
    X(datos.Ncol) = X(datos.Ncol-1) + dx(datos.Ncol)/2;
    
    Y (1) = 0; Y(2) = dy(1)/2;
    for j = 3:datos.Nrow-1
        Y(j) = Y(j-1) + dy(j-1);
    end
    Y(datos.Nrow) = Y(datos.Nrow-1) + dy(datos.Nrow)/2;
    
    for i=1:datos.Nrow
        for j=1:datos.Ncol
            Coll_X(i,j) = X(i);
            Coll_Y(i,j) = Y(j);
        end
    end

    
function [ X, Y ] = STAGG_MESH_X ( datos, C ) 
    
    vector_x = zeros( datos.Nx+1 ,1 );
    
    vector_x(1) = 0; % La malla comença a la paret
    
    for i = 2:datos.Nx
        vector_x(i) = vector_x(i-1) + C.dx(i);
    end
    
    vector_x(datos.Nx+1) = datos.L;
    
    for i = 1:datos.Nrow_sx+1
       for j = 1:datos.Ncol_sx
            X(i,j) = vector_x(j);
            Y(i,j) = C.Y(i);
       end
    end
    
    fprintf('    Staggered Mesh X calculada\n');

function [ X, Y ] = STAGG_MESH_Y ( datos, C ) 
    
    vector_y = zeros( datos.Ncol_sy ,1 );
    
    vector_y(1) = 0; % La malla comença a la paret
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dx(i);
    end
    
    vector_y(datos.Ny+1) = datos.L;
    
    for i = 1:datos.Nrow_sy
       for j = 1:datos.Ncol_sy+1
            X(i,j) = C.X(j);
            Y(i,j) = vector_y(i);
       end
    end
    
    fprintf('    Staggered Mesh Y calculada\n');


function delta = spacing_no_uniform (gamma, L, i, N, uniform )

    if uniform == true
        delta = L * i / N;
    else
        delta = L*0.5*(1.0+tanh(gamma*(2.0*i/N-1.0))/tanh(gamma));
    end
    
function S = CondicionesIniciales (datos, C)
    
    fprintf('\n-( 2 )- Condiciones iniciales establecidas \n');
    
    S.P = zeros(datos.Nrow, datos.Ncol);
    S.U = zeros(datos.Nrow_sx, datos.Ncol_sx);
    S.V = zeros(datos.Nrow_sy, datos.Ncol_sy);
    
    S.P_prev = ones(datos.Nrow, datos.Ncol);
    S.U_prev = zeros(datos.Nrow_sx, datos.Ncol_sx);
    S.V_prev = zeros(datos.Nrow_sy, datos.Ncol_sy);
    
    S.R_u = zeros(datos.Nrow_sx, datos.Ncol_sx);
    S.R_v = zeros(datos.Nrow_sy, datos.Ncol_sy);
    S.R_u_prev = zeros(datos.Nrow_sx, datos.Ncol_sx);
    S.R_v_prev = zeros(datos.Nrow_sy, datos.Ncol_sy);
    
    S.time = 0;
    S.dt = datos.dt_inicial;

    
 function S = Solve_NS (datos, C);
    
    fprintf('\n-( 3 )- Calculando solución... \n');
     
    converg = true; % no vull que entri encara
    
    S = 0;
    %{
    while converg == false || S.time < datos.tfinal 
         
        % Update de variables
            S.U_prev = S.U;
            S.V_prev = S.V;
        
        % Factor R
            [S.R_u s.R_v] = Compute_R(datos,C,S);
   
        
    end
    %}
    
    % Print results
function [R_u R_v] = Compute_R (datos,C,S)

    
    % Update de variables
        S.R_u_prev = S.R_u;
        S.R_v_prev = S.R_v;
        
   %% NOMENCLATURA
    % Fe = flux sortint massic en la cara e (EAST)
    % uE = velocitat direccio x en el node E (EAST)
    % ue = velocitat direccio x en la cara e (EAST)
    % Se = superficie de la cara e (EAST)
    
    % Vegeu que sub index e indica cara; E indica node
    
    for i=1:datos.Nrow_sx
        for j=1:datos.Ncol_sx
            
            UE = S.U(i,j+1);
            UW = S.U(i,j-1);
            UN = S.U(i+1,j);
            US = S.U(i-1,j);
            UP = S.U(i,j+1);
            
            Se = C.stagX_y(i,j+1) - C.stagX_y(i,j);
            Sw = C.stagX_y(i,j) - C.stagX_y(i,j);
            Sn = C.stagX_x(i,j) - C.stagX_x(i,j);
            Ss = C.stagX_x(i,j) - C.stagX_x(i,j);
            
            % Aproximacio terme convectiu
            ue = (U(i,j+1)+U(i,j))/2; Fe = ue*Se;
            uw = (U(i,j+1)+U(i,j))/2; Fw = uw*Sw
            un = (U(i,j+1)+U(i,j))/2;
            us = (U(i,j+1)+U(i,j))/2;
            vn = (V(i,j+1)+V(i,j))/2; Fn = vn*Sn;
            vs = (V(i,j+1)+V(i,j))/2; Fs = vs*Ss;
            
        end
    end
        
         
       
        