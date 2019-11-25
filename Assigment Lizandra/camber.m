%% Function Camber %%
% calculo ordenada Z NACA 4 digitos
% perfil cuerda unitaria!

function z = camber (x,p,m)
    if x<p
        z=(m/p^2)*(2*p*x-x^2);
    else
        z=(m/(1-p)^2)*(1-2*p+2*p*x-x^2);
    end
end
