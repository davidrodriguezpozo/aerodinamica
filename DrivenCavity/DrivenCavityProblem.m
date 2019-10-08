function main

clear all
clc 
close all

Vx = 10; Vx = 10;
datos = Driven_Input(Vx,Vx);
C = Driven_mesh (datos); 

disp('A calcular')



function  C = Driven_mesh (datos)
    
    fprintf('-( 1 )- Calculando mallas...\n'); 
        
    % Preallocation of variables
    
    C.X = zeros (datos.Nx,datos.Ny);
    C.Y = zeros (datos.Nx,datos.Ny);    
    C.dx = zeros (datos.Nx,1);
    C.dy = zeros (datos.Ny,1);
    
    C.stagX_x = zeros (datos.Nx,datos.Ny);
    C.stagX_y = zeros (datos.Nx,datos.Ny);
    
    C.stagY_x = zeros (datos.Nx,datos.Ny);
    C.stagY_y = zeros (datos.Nx,datos.Ny);
    
    [ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = Driven_COLLOCATED_MESH ( datos );

    [ C.stagX_x, C.stagX_y ] = Driven_STAGG_MESH_X ( datos, C );
    
    [ C.stagY_x, C.stagY_y ] = Driven_STAGG_MESH_Y ( datos, C );

    
    
function [ X, Y, dx, dy, Coll_X, Coll_Y ] = Driven_COLLOCATED_MESH ( datos )
    
    % Delta X
    for i = 1:datos.Nx
            dx(i) = datos.L * (i+1) / datos.Nx - datos.L * i / datos.Nx ;
    end 
    
    % Delta Y
    for j = 1:datos.Ny
            dy(j) = datos.H * (j+1) / datos.Ny - datos.H * j / datos.Ny ;
    end 
    
    % COLLOCATED MESHES
    
    X(1) = dx(1)/2;
    for i = 2:datos.Nx
        X(i) = X(i-1) + dx(i-1);
    end
    
    Y(1) = dy(1)/2;
    for j = 2:datos.Ny
        Y(j) = Y(j-1) + dy(j-1);
    end
    
    for i=1:datos.Nx
        for j=1:datos.Ny
            Coll_X(i,j) = X(j);
            Coll_Y(i,j) = Y(i);
        end
    end
    
    
function [ X, Y ] = Driven_STAGG_MESH_X ( datos, C ) 
    
    vector_x = zeros( datos.Nx ,1 );
    
    vector_x(1) = 0; % La malla comen�a a la paret
    
    for i = 2:datos.Nx
        vector_x(i) = vector_x(i-1) + C.dx(i);
    end
    
    vector_x(datos.Nx+1) = datos.L;
    
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = vector_x(i);
            Y(i,j) = C.Y(j);
       end
    end
    
    fprintf('    Staggered Mesh X calculada\n');
    
    
    function [ X, Y ] = Driven_STAGG_MESH_Y ( datos, C ) 
    
    vector_y = zeros( datos.Ny ,1 );
    
    vector_y(1) = 0; % La malla comen�a a la paret
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dy(i);
    end
    
    vector_y(datos.Ny) = datos.L;
    
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = C.X(i);
            Y(i,j) = vector_y(j);
       end
    end
 
    fprintf('    Staggered Mesh Y calculada\n');
 
function datos = Driven_Input(Vx,Vy)

    datos.Nx = Vx+2;
    datos.Ny = Vy+2;
    
    datos.Vx = Vx;
    datos.Vy = Vy;
    
    datos.L = 1;
    datos.H = 1;


 
