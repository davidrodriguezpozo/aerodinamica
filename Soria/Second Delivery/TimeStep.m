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

function delta_t = TimeStep (datos,u,v,delta_t)
%     %Time for convective term
%     delta = datos.L/datos.Vx;
%     Tc_u = min((min(delta./abs(u))));
%     Tc_v = min((min(delta./abs(v))));
%     Tc = min (Tc_u,Tc_v);
%     %Time for diffusive term
%     nu = datos.mu/datos.rho;
%     Td = 0.5*(delta^2/nu);
%     
%     %Step time:
%     delta_T = 0.001 * min(Tc,Td);
%     %delta_T =0.01 ;
volumes_y = datos.Vy;
volumes_x = datos.Vx;
delta = datos.L/datos.Vx;
delta_tc_u_min=100;
delta_tc_v_min=100;
delta_td = 0.5*delta^2*datos.rho/datos.mu;

for i=2:volumes_y+1
    for j=2:volumes_x+1
       
     if (abs((delta)/(u(i,j)))<delta_tc_u_min)
         
         delta_tc_u_min=abs((delta)/(u(i,j)));
         
     end
     
      if (abs((delta)/(v(i,j)))<delta_tc_v_min)
         
         delta_tc_v_min=abs((delta)/(v(i,j)));
         
      end
     
    end    
end
     
delta_tc = min([delta_tc_v_min, delta_tc_u_min]);

if (0.2*min([delta_tc, delta_td])<delta_t)

    delta_t=0.2*min([delta_tc, delta_td]);
    
end
    delta_t=0.2*min([delta_tc, delta_td]);
%delta_t = 0.01;
end
