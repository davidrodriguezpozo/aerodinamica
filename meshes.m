% MESHES

function meshes

close all 
clc 
clear all

    datos.Vx = 3;
    datos.Vy = 3;
    datos.Nx = datos.Vx+2;
    datos.Ny = datos.Vy+2;
    
    datos.L = 1;
    datos.H = 1;
    
    datos.uniform = true;
    datos.gamma = 1;
    
    malla = 3;
    
    C = MESH(datos,malla);
      
    
function    C = MESH (datos, malla)
    
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
    
    [ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( datos );

    [ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( datos, C );
    
    [ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( datos, C );
    
    figure; hold on;
    plot(C.Coll_x , C.Coll_y,'k.','MarkerSize', 35) ; 
    plot(C.stagX_x , C.stagX_y,'b.','MarkerSize',20) ;
    plot(C.stagY_x , C.stagY_y,'g.') ;
    ylim ([-0.1 1.1])
    xlim ([-0.1 1.1])
    
    % Collocated Mesh

    V1x = [C.dx(1) datos.L-C.dx(datos.Nx)];
    V1y = [C.dy(1) datos.H-C.dy(datos.Ny) ]; 
    V2x = [ datos.L-C.dx(datos.Nx)   datos.L-C.dx(datos.Nx) ]; 
    V2y = [ datos.H-C.dy(datos.Ny)   datos.H-C.dy(datos.Ny)];
    V3x = [C.dx(1) C.dx(1)];
    V3y = [C.dy(1) C.dy(1)];
   
    V5 = [C.dx(1)+C.dx(2) C.dx(1)+C.dx(2)];
    V6 = [C.dy(1) datos.H-C.dy(datos.Ny) ]; 
    V7 = [C.dx(1)+C.dx(2)+C.dx(3) C.dx(1)+C.dx(2)+C.dx(3)];
    
    V8 = [C.dx(1) datos.H-C.dy(datos.Ny) ];
    V9 = [C.dy(1)+C.dy(2) C.dy(1)+C.dy(2) ]; 
    V10 = [C.dy(1)+C.dy(2)+C.dy(3) C.dy(1)+C.dy(2)+C.dy(3)];
    
    if malla == 1
        hold on;
        plot(V1x,V2y,'k');  plot(V3x,V1y,'k'); plot(V2x,V1y,'k'); plot(V1x,V3y,'k');
        plot(V5,V6,'-.k'); plot(V7,V6,'-.k');
        plot(V8,V9,'-.k'); plot(V8,V10,'-.k');
        title('Collocated mesh'); 
    end
    
    % Staggered Mesh X
    if malla == 2
        V1x = V1x - C.dx(1)/2; V3x = V3x - C.dx(1)/2;
        V2x = V2x - C.dx(1)/2; 
    
        plot(V1x,V2y,'b');  plot(V3x,V1y,'b'); plot(V2x,V1y,'b'); plot(V1x,V3y,'b');
        plot(V5- C.dx(1)/2,V6,'-.b'); plot(V7- C.dx(1)/2,V6,'-.b');
        plot(V8- C.dx(1)/2,V9,'-.b'); plot(V8- C.dx(1)/2,V10,'-.b');
        title('Staggered mesh x-direction'); 
    end
        
    % Staggered Mesh Y
    if malla == 3
        V1y = V1y - C.dy(1)/2; V3y = V3y - C.dy(1)/2;
        V2y = V2y - C.dy(1)/2; 
        
        plot(V1x,V2y,'g');  plot(V3x,V1y,'g'); plot(V2x,V1y,'g'); plot(V1x,V3y,'g');
        plot(V5,V6- C.dy(1)/2,'-.g'); plot(V7,V6- C.dy(1)/2,'-.g');
        plot(V8,V9- C.dy(1)/2,'-.g'); plot(V8,V10- C.dy(1)/2,'-.g');
        title('Staggered mesh y-direction'); 
    end
    
    fprintf('---> Mallas calculada\n')      
   
    
function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( datos )
    
    % Delta X
    for i = 1:datos.Nx
            dx(i) = spacing_no_uniform( datos.gamma, datos.L, i, datos.Nx, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.L, i-1, datos.Nx, datos.uniform );
    end 
    
    % Delta Y
    for j = 1:datos.Ny
            dy(j) = spacing_no_uniform( datos.gamma, datos.H, j, datos.Ny, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.H, j-1, datos.Ny, datos.uniform );
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
    
    
function [ X, Y ] = STAGG_MESH_X ( datos, C ) 
    
    vector_x = zeros( datos.Nx ,1 );
    
    vector_x(1) = 0; % La malla comença a la paret
    
    for i = 2:datos.Nx
        vector_x(i) = vector_x(i-1) + C.dx(i);
    end
    
    %vector_x(datos.Nx+1) = datos.L;
    
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = vector_x(i);
            Y(i,j) = C.Y(j);
       end
    end
    
    fprintf('    Staggered Mesh X calculada\n');
    
    
    function [ X, Y ] = STAGG_MESH_Y ( datos, C ) 
    
    vector_y = zeros( datos.Ny ,1 );
    
    vector_y(1) = 0; % La malla comença a la paret
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dy(i);
    end
    
    %vector_y(datos.Ny) = datos.L;
    
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = C.X(i);
            Y(i,j) = vector_y(j);
       end
    end
    
    fprintf('    Staggered Mesh Y calculada\n');

    
function delta = spacing_no_uniform (gamma, L, i, N, uniform )

    if uniform == true
        delta = L * i / N;
    else
        delta = L*0.5*(1.0+tanh(gamma*(2.0*i/N-1.0))/tanh(gamma));
    end
    