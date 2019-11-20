function main

clear all
clc 
close all

N = 10; 
Re = 100;

datos = Driven_Input(N,Re);
C = Driven_mesh (datos); 

3disp('A calcular')



function  C = Driven_mesh (datos)
    
    fprintf('-( 1 )- Calculando mallas...\n'); 
        
    % Preallocation of variables
    
    C.X = zeros (datos.N,datos.N);
    C.Y = zeros (datos.N,datos.N);    
    
    C.stagX_x = zeros (datos.Nx_sx,datos.Ny_sx);
    C.stagX_y = zeros (datos.Nx_sx,datos.Ny_sx);
    
    C.stagY_x = zeros (datos.Nx_sy,datos.Ny_sy);
    C.stagY_y = zeros (datos.Nx_sy,datos.Ny_sy);
    
    [ C.X, C.Y ] = Driven_COLLOCATED_MESH ( datos );

    [ C.stagX_x, C.stagX_y ] = Driven_STAGG_MESH_X ( datos, C );
    
    [ C.stagY_x, C.stagY_y ] = Driven_STAGG_MESH_Y ( datos, C );

    
    
function [ X, Y ] = Driven_COLLOCATED_MESH ( datos )
    
    N = datos.N;
    dx = datos.L/N; dy = datos.H/N;
    
    X(1)=0; X(2) = dx/2; 
    for i=3:datos.Nx-1
        X(i) = X(i-1) + dx; 
    end
    X(datos.Nx) = X(datos.Nx-1) + dx/2; 
    
    Y(1)=0; Y(2) = dy/2; 
    for j=3:datos.Ny-1
        Y(j) = Y(j-1) + dy; 
    end
    Y(datos.Ny) = Y(datos.Ny-1) + dy/2; 
    
    fprintf('    Collocated Mesh calculada\n');
    
    
    
function [ X, Y ] = Driven_STAGG_MESH_X ( datos, C ) 
    
    N = datos.N;
    dx = datos.L/datos.N;
    
    vector_x = zeros(datos.N,1);
    
    vector_x(1) = 0; % La malla comen�a a la paret
    
    for i = 2:datos.N
        vector_x(i) = vector_x(i-1) + dx;
    end
    
    vector_x(N+1) = datos.L;
    
    for i = 1:datos.Nx_sx
       for j = 1:datos.Ny_sx
            X(i,j) = vector_x(i);
            Y(i,j) = C.Y(j);
       end
    end
    
    fprintf('    Staggered Mesh X calculada\n');
    
    
function [ X, Y ] = Driven_STAGG_MESH_Y ( datos, C ) 
    
    N = datos.N;
    dy = datos.H/datos.N;
    
    vector_y = zeros(datos.N,1);
    
    vector_y(1) = 0; % La malla comen�a a la paret
    
    for j = 2:datos.N
        vector_y(j) = vector_y(j-1) + dy;
    end
    
    vector_y(N+1) = datos.H;
    
    for i = 1:datos.Nx_sy
       for j = 1:datos.Ny_sy
            X(i,j) = C.X(i);
            Y(i,j) = vector_y(j);
       end
    end
    
    fprintf('    Staggered Mesh Y calculada\n');
    
    
function datos = Driven_Input(N,Re)
  
    datos.N = N;
    
    datos.Nx = N+2;
    datos.Ny = N+2;
    
    datos.Nx_sx = datos.Nx -2;
    datos.Ny_sx = datos.Ny -1;
    
    datos.Nx_sy = datos.Nx -1;
    datos.Ny_sy = datos.Ny -2;
    
    datos.L = 1;
    datos.H = 1;
    datos.rho = 1;
    datos.delta = 10e-6;
    datos.delta2 = 10e-6;
    datos.Re = Re;
    datos.final_time = 10;
    
    
    


 
