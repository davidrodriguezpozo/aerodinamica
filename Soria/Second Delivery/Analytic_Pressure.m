%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    FUNCTION ANALYTIC PRESSURE     %%%%%%%%%%%%%
%%% This function computes the analytical solution for  %%%
%%% the pressure defined in this own function. The solu-%%%
%%% tion is obtained for the x and y directions.        %%%
%%% Initial inputs are used as a initial data.          %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function p = Analytic_Pressure (datos, C, time)

    for i = 1:datos.Nx
        for j=1:datos.Ny
        
        mu = datos.mu;
        rho = datos.rho;
        x = C.X(i);    
        y = C.Y(j);
        p(i,j) = -(exp(-8*pi^2*mu*time))^2 * rho * (cos(4*pi*x)+ cos(4*pi*y))/4;
           
        end 
    end
