function main_B

clear all
clc
close all

%Es la mateixa funcio que el main pero sense les divisions de VC pel M.M.S.

%{  
            El Codi hauria de fer aix� 

 1. Timestep for stability
 2. u_n_1 i u_n conegudes --> R_n_1 i R_n
 3. Evaluar u_p = u_n + ?t(3/2 * R_n - 1/2 * R_N_1
        R = - Conv(u)/V + Diff(u)/V
 4. Resoldre l'equaci� de Poisson A * pseudo_p = Div(u_p)
 5. Obtenir pseudo_p
 6. Obtenir u_n1 --> u_n1 = u_p - grad(pseudo_p) 
 7. timestep + 1 
 %}

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 


Vx = 10; %N� de divisions en x de V.C.
Vy = 10;

datos = Input(Vx,Vy);
C = meshes (datos, datos.malla);

matriu_A = A_laplace(datos,C); %Pseudo-pressure matrix A.

v = zeros(Vx,Vy);
u = zeros(Vx,Vy);
u(2,2) = 1;
[nodal_mesh num] = nodalmesh(Vx,Vy)

u_p = divergencia_u(datos, u, v, nodal_mesh, num);

pseudo_p = zeros(Vx*Vx,1); %Pseudo-pressure vector NOT KNOWN 

%Per solucionar el problema de la matriu A singular (no es podria invertir)

matriu_A(1,1)=-5;

pseudo_p = inv(matriu_A)*u_p;

[p_gradX p_gradY] = gradient_p(datos, pseudo_p, nodal_mesh, num)

u_n1 = u - p_gradX;
v_n1 = v - p_gradY;

divvvvvv = divergencia_u(datos, u_n1, v_n1, nodal_mesh, num);

%UEEEEEEA dona zero la div jejeje
%  |
%  |
%  V
%  V
sum(divvvvvv)

%% AÑADO
    rho = 1;
    
    delta_T = 1;
    
    p = pseudo_p * delta_T / rho; 
    
    k = 1; t = 1;
    for i=1:size(p,1)
        if k <= datos.Vx
            p_mat(k,t) = p(i);
            k = k+1;
        else
            k = 1; t = t+1;
            p_mat(k,t) = p(i);
        end
    end
    
    disp('finished');



