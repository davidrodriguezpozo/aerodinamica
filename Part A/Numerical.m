  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         NUMERICAL FUNCTION        %%%%%%%%%%%%%
%%% This function computes the Numerical solution for   %%%
%%% the convective and the diffusive term for the x-di- %%%
%%% rection and the y-direction.                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v)
% The inputs: the initial values (datos), the meshes (C), and the velocity field computed in other function. 
% The outputs: the convectives and diffusives terms. 

        conv_u = convection(u,v,datos.L,datos.H);
        diff_u = diffusion_u(u,datos.L);

        conv_v = convection(v,u,datos.H,datos.L);
        diff_v = diffusion_u(v,datos.H);
