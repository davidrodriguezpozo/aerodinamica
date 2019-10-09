function main_A

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 

div = [3 5 10 20 30 40 50 60 70 80 90 100];

for k = 1:length(div)
    
Nx_vector(k) = div(k); 
Ny_vector(k) = div(k);

Vx = Nx_vector(k);

datos = Input(Vx,Vx);
C = meshes (datos, datos.malla);

% FEM SOLUCIO ANALITICA
malla = 'x';
[S.cu_anal S.du_anal] = Analytic (datos, C, malla);

malla = 'y';
[S.cv_anal S.dv_anal] = Analytic (datos, C, malla);

%disp('Analitica calculada');

clear malla

% PERFIL VELOCITATS
[u v ] = Velocitats (datos, C);

[S.cu S.diffu S.cv S.diffv] = Numerical (datos, C, u, v);



if k == length(div)   
    plot_Perfil_Vel (datos,C,u,v)   
end

   [error_du(k) error_cu(k) error_dv(k) error_cv(k)] = Compute_error ( S, k);

end

plotsA(Nx_vector, error_du, error_dv, error_cu, error_cv);







 

    

   
   