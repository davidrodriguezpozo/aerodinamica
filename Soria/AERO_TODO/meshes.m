%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION MESHES         %%%%%%%%%%%%%
%%% This function allows the programer to allocate the  %%%
%%% different points of the three meshes needed to sol- %%%
%%% ve the NS equations. Furthermore, three functions   %%%
%%% are contempled in order to fullfill all the require-%%%
%%% ments. First of all, a Colloated-Mesh function will %%% 
%%% be performed moving the y-axis center and the x-axis%%%
%%% center up and right, respectively. Afterwards, two  %%%
%%% Staggered-Mesh functions will be executed (from x & %%%
%%% y-axis) building two other meshes into a "regular"  %%%
%%% performing position (they are explained below). Fin-%%%
%%% ally a position-incremental function is defined and %%%
%%% the sum of all is explained in the meshes function, %%%
%%% that works as a "main".                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function  C = meshes (datos, malla)
    % A C class is defined in order to organise all the outputs of the
    % function.
   
        
    % Preallocation of variables: X and Y positions (corresponding to
    % collocated locations), incrementals, staggered (x and y) matrixes
    % which will allocate the nodal positions.
    C.X = zeros (datos.Nx,datos.Ny);
    C.Y = zeros (datos.Nx,datos.Ny);    
    C.dx = zeros (datos.Nx,1);
    C.dy = zeros (datos.Ny,1);
    
    % Inicializtion of:
    C.stagX_x = zeros (datos.Nx,datos.Ny); %X-staggered mesh for x-axis
    C.stagX_y = zeros (datos.Nx,datos.Ny); %X-staggered mesh for y-axis
    C.stagY_x = zeros (datos.Nx,datos.Ny); %Y-staggered mesh for x-axis
    C.stagY_y = zeros (datos.Nx,datos.Ny); %Y-staggered mesh for y-axis
    
    %Filling the meshes with each function
    [ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( datos );
    [ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( datos, C );
    [ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( datos, C );
    
    
    
  

% Collocated_mesh function build the collocated mesh with no uniform meshes
% (computing delta X and delta Y for every type of distribution)
function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( datos )
    
    dx = datos.L/datos.Vx;
    
    dy = datos.H/datos.Vy;
    
    % Computing the X-position of COLLOCATED MESH
    X(1) = -dx/2; % X is inicialized in the average center of both VC faces
    for i = 2:datos.Nx
        X(i) = X(i-1) + dx; % Delta X is added in each defined position from the second one
    end
    
    
    % Computing the Y-position of COLLOCATED MESH
    Y(1) = -dy/2; % Y is inicialized in the average center of both VC faces
    for j = 2:datos.Ny
        Y(j) = Y(j-1) + dy; % Delta Y is added in each defined position from the second one
    end
    
    % Finally, the whole mesh is joined using as reference the coordenate
    % axis. 
    for i=1:datos.Nx
        for j=1:datos.Ny
            Coll_X(i,j) = X(j);
            Coll_Y(i,j) = Y(i);
        end
    end


% Stagg_mesh_x function build the staggered-x mesh with no uniform meshes
% (computing delta X for every type of distribution)
function [ X, Y ] = STAGG_MESH_X ( datos, C ) 
    
    % A x-direction vector is inicialized (with Nx size)
    vector_x = zeros( datos.Nx ,1 );
    
    % The mesh starts in the first position (0)
    vector_x(1) = 0;
    
    for i = 2:datos.Nx
        vector_x(i) = vector_x(i-1) + C.dx; 
        % From the first position, the incremental calculated before is 
        % added in each point. Adding the same delta that in the collocated
        % mesh assures the same VC length.
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %vector_x(datos.Nx+1) = datos.L;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % The X-Position is computed related to a vertical position. 
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = vector_x(j);
            Y(i,j) = C.Y(i);
       end
    end
   
    
    
    function [ X, Y ] = STAGG_MESH_Y ( datos, C )
    % This function works as the function STAGG_MESH_X
    
    vector_y = zeros( datos.Ny ,1 );
    
    vector_y(1) = 0; 
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dy;
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %vector_y(datos.Ny) = datos.L;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = C.X(j);
            Y(i,j) = vector_y(i);
       end
    end
