%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION ERROR          %%%%%%%%%%%%%
%%% Comparing the analytic and the numeric results of   %%%
%%% the convective and diffusive terms, this function   %%%
%%% obtains the relative error between each measure,    %%%
%%% allowing the designer to plot the variation among   %%%
%%% an incremental variable (Mesh size or time, p.e.).  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error = ERROR(numeric,analytic)
    % Firstly the error is inisialised. 
    error = 0;
    % For each point of the domain, the error is calculated. 
    for i = 1:size(numeric,1)
       for j = 1:size(numeric,2)
          
               a = abs(analytic(i,j)-numeric(i,j));
         
            if a > error 
                error = a;
            end
       end
    end
