%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      FUNCTION DELTA_T               %%%%%%%%%%%
%%% This function computes the time-step incremental for%%%
%%% the part-C soltion dividing into the convective and %%%
%%% diffusive time-step. The minimum between both is    %%%
%%% used in other functions.                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function delta_T = TimeStep (datos,u,v)
        %Time for convective term
       delta = datos.L/datos.Vx;
    Tc_u = min(datos.L/(datos.Vx*max(max(u))));
    Tc_v = min(datos.H/(datos.Vy*max(max(v))));
    Tc = min (Tc_u,Tc_v);
    %Time for diffusive term
    Td = 0.5*(datos.L/datos.Vx)*(datos.H/datos.Vy) / datos.mu;
    
    %Step time:
    delta_T = 0.2 * min(Tc,Td);
    % delta_T =0.02 ;
