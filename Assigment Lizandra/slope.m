%% Function Slope %%
% calculo de la pendiente NACA 4 digitos

function dzdx = slope (x,p,m)
    if x<p
        dzdx=(2*m/p^2)*(p-x);
    else
        dzdx=(2*m/(1-p)^2)*(p-x);
    end
end
