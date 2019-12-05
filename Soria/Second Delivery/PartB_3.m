

function PartB_3

    clear all 
    clc
    close all 

%     Vx = [3 5 10 20 30 40]; 
%     Vy = [3 5 10 20 30 40];
    Vx = 3;
    Vy = 3;
    Re = 1000;
    
%     for i=1:length(Vx)
%         datos = INPUT(Vx(i),Vy(i),Re);
        datos = INPUT(Vx,Vy,Re);
        datos.F = 1;
        
        %Random predicted velocity field that satisfies mass conservation
        u_p = zeros(datos.Nx,datos.Ny); 
        v_p = zeros(datos.Nx,datos.Ny);
        u_p(2,2) = 5;
        u_p = haloupdate(u_p);
        v_p = haloupdate(v_p);
        C = meshes (datos, datos.malla);  
        [u, v, P, time] = SolverB(datos, C, u_p, v_p);
        display(u);
        display(v);
        latex(sym(u_p))
        latex(sym(vpa(u)))
        latex(sym(vpa(v)))
%     end
    
    %PlotsB (C,P,Vx,errorP)
    
    
    
    
       

