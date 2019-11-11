

function main_B

    clear all 
    clc
    close all 

    Vx = [3 5 10 20 30 40]; 
    Vy = [3 5 10 20 30 40]; 
    Re = 1000;
    
    for i=1:length(Vx)
    
    datos = INPUT(Vx(i),Vy(i),Re);
    datos.F = 1;
    C = meshes (datos, datos.malla);  
    [u, v, P, time, errorP(i)] = SolverB(datos, C);
    
    end
    
    PlotsB (C,P,Vx,errorP)
    
    
    
    
       

