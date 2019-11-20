  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         FUNCTION VELOCITATS         %%%%%%%%%%%
%%% This function, obviously, iterates between the stag-%%%
%%% gered meshes computing the velocities in the face of%%%
%%% this kind of meshes. We define the meshes using x or%%%
%%% y values, imposing the condition in the main func-  %%%
%%% tion.                                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u v ] = Velocitats (datos, C)
% The inputs: initial values (datos) and meshes definitons (c).
% The outputs: u is the velocity for x-axis and v for y-axis.

% Starting the double loop iteration for the x and y meshes-dimensions. 
for i = 1:datos.Nx
    for j = 1:datos.Ny
        % For x or y staggered meshes. 
        x = C.stagX_x(i,j);
        y = C.stagX_y(i,j);
        % The velocity matrix in the x direction and staggered mesh is refilled.
        u(i,j) = datos.F*cos(2*pi*x)*sin(2*pi*y);
        
        % Item for the part above. 
        x = C.stagY_x(i,j);
        y = C.stagY_y(i,j);
        
        v(i,j) = -datos.F*cos(2*pi*y)*sin(2*pi*x);
        
    end
end
