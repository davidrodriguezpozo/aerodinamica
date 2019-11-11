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
    
    % Computing Delta X
    for i = 1:datos.Nx
            dx(i) = spacing_no_uniform( datos.gamma, datos.L, i, datos.Nx, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.L, i-1, datos.Nx, datos.uniform );
    end 
    
    % Computing Delta Y
    for j = 1:datos.Ny
            dy(j) = spacing_no_uniform( datos.gamma, datos.H, j, datos.Ny, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.H, j-1, datos.Ny, datos.uniform );
    end 
    
    % Computing the X-position of COLLOCATED MESH
    X(1) = dx(1)/2; % X is inicialized in the average center of both VC faces
    for i = 2:datos.Nx
        X(i) = X(i-1) + dx(i-1); % Delta X is added in each defined position from the second one
    end
    
    
    % Computing the Y-position of COLLOCATED MESH
    Y(1) = dy(1)/2; % Y is inicialized in the average center of both VC faces
    for j = 2:datos.Ny
        Y(j) = Y(j-1) + dy(j-1); % Delta Y is added in each defined position from the second one
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
        vector_x(i) = vector_x(i-1) + C.dx(i); 
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
            X(i,j) = vector_x(i);
            Y(i,j) = C.Y(j);
       end
    end
   
    
    
    function [ X, Y ] = STAGG_MESH_Y ( datos, C )
    % This function works as the function STAGG_MESH_X
    
    vector_y = zeros( datos.Ny ,1 );
    
    vector_y(1) = 0; 
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dy(i);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %vector_y(datos.Ny) = datos.L;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = C.X(i);
            Y(i,j) = vector_y(j);
       end
    end
    
   

% Spaging_no_uniform computes the Deltas (separations between VC faces)
% for every type of distribution    
function delta = spacing_no_uniform (gamma, L, i, N, uniform )

    % The reader should note that this function needs to know if the mesh
    % is uniform or not, creating a more efficient solution. 
    if uniform == true
        % L is the total mesh length. i is the position. N is the number of
        % nodes.
        delta = L * i / N;
    else
        delta = L*0.5*(1.0+tanh(gamma*(2.0*i/N-1.0))/tanh(gamma));
    end
    
