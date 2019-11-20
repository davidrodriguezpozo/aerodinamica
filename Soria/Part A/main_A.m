%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION MAIN           %%%%%%%%%%%%%
%%% This function sumarises all the the functions deve- %%%
%%% loped in the code. As it will be showed, it follows %%%
%%% a logic path in order to solve every step in the    %%%
%%% initial problem. First of all, a mesh size is defi- %%%
%%% ned (in this case, a growing loop is used). Initia- %%%
%%% lizing the input values, it is possible to define a %%%
%%% mesh for every part of the problem. Then, an analyt-%%%
%%% ical solution and a numerical solution are computed %%%
%%% preparing the final step, which is reporting the    %%%
%%% error between both solutions. Finally, a plot of    %%%
%%% every iteration error is done showing the optimal   %%%
%%% solution (because every iteration implies differents%%%
%%% sizes of mesh.                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function main_A

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 

% Starting of the iteration
div = [3 5 10 20 30 40 50 60 70 80 90 100];

for k = 1:length(div)
% Defining the meshes dimensions    
Nx_vector(k) = div(k); 
Ny_vector(k) = div(k);

Vx = Nx_vector(k);

% Defining the "datos" structure, where the inputs are kept.
datos = Input(Vx,Vx);
C = meshes (datos, datos.malla);

% Initializing functions which will calculate the analyitical solution
malla = 'x'; % It will allow to separate problems between both axis. 
[S.cu_analytic S.du_analytic] = Analytic (datos, C, malla);

malla = 'y';
[S.cv_analytic S.dv_analytic] = Analytic (datos, C, malla);

clear malla

% This function computes the velocity field in both directions with the initial values. 
[u v ] = Velocitats (datos, C);

[S.cu S.diffu S.cv S.diffv] = Numerical (datos, C, u, v);

% Then, we show the velocity field.
if k == length(div)   
    plot_Perfil_Vel (datos,C,u,v)   
end

% Finally, the error is computed for the diffusive and the convective terms for both directions in every iteration.
   [error_du(k) error_cu(k) error_dv(k) error_cv(k)] = Compute_error ( S, k);

end

% For better understanding, the error is computed depending on the mesh dimension. 
plotsA(Nx_vector, error_du, error_dv, error_cu, error_cv);







 

    

   
   
