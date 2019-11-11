%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION NUMERICAL       %%%%%%%%%%%%
%%% This function just executes the functions which so- %%%
%%% lute the convection and diffusive term.             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v)

        [conv_u conv_v] = convection_main(u,v,datos.L,datos.H,C);
        diff_u = diffusion_u(u,datos.L);
        diff_v = diffusion_u(v,datos.H);
