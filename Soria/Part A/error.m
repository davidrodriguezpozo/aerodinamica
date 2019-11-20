  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION ERROR          %%%%%%%%%%%%%
%%% This function analyse the error between the numeric-%%%
%%% al and the analytical solutions of the convective   %%%
%%% and the diffusive terms of the NS equation. This job%%%
%%% will be done substrating analitical and numerical   %%%
%%% terms of each solution.                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function error = ERROR(numeric,analytic)
 % The inputs: numeric and analytic terms of diffusive and convective solutions.
 % The outputs: the error.
 
 % Firstly, the error is initialized in case of avoid weird results.
    error = 0;

 % In this double bucle, the developer intends to analyse the error in each position of the mesh. 
    for i = 1:size(numeric,1)
       for j = 1:size(numeric,2)
               a = abs(analytic(i,j)-numeric(i,j));
            if a > error 
                error = a;
            end
       end
    end
