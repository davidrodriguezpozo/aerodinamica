function main_B_2

%Es la mateixa funcio que el main pero sense les divisions de VC pel M.M.S.

%{  
                                El Codi hauria de fer aixï¿½ 

 1. Timestep for stability
 2. u_n_1 i u_n conegudes --> R_n_1 i R_n
 3. Evaluar u_p = u_n + ?t(3/2 * R_n - 1/2 * R_N_1
        R = - Conv(u)/V + Diff(u)/V
 4. Resoldre l'equaciï¿½ de Poisson A * pseudo_p = Div(u_p)
 5. Obtenir pseudo_p
 6. Obtenir u_n1 --> u_n1 = u_p - grad(pseudo_p) 
 7. timestep + 1 
 %}

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

clear all 
clc
close all 



Vx = 3; %nombre de divisions en x de V.C SENSE HALO
Vy = 3;
R_uant=zeros(Vx+2, Vx+2);
R_vant=zeros(Vx+2, Vx+2);
R_u = zeros(Vx+2, Vx+2);
R_v = zeros(Vx+2, Vx+2);
Re = 10000;
datos = INPUT(Vx,Vy,Re);
C = meshes (datos, datos.malla);

v = zeros(datos.Nx,datos.Ny);
u = zeros(datos.Nx,datos.Ny);
u(3,3) = 1; % El camp de velocitats té HALO
u = haloupdate(u);
%[nodal_mesh num] = nodalmesh(Vx,Vy)
[nodal_mesh, num] = nodalmesh(Vx,Vy);
%delta_t = TimeStep(u,v);
delta_t = 0.02;

[conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);

[u_p, v_p, R_uant, R_vant] = predicted(u,v,delta_t,datos, conv_u, diff_u, conv_v, diff_v, R_uant, R_vant,R_u,R_v);

div_u_p = divergencia_u(datos, u_p, v_p, nodal_mesh, num); %u_p = [Vx*Vy,1]

pseudo_p = zeros(Vx*Vx,1); %Pseudo-pressure vector NOT KNOWN 


%Per solucionar el problema de la matriu A singular (no es podria invertir)

matriu_A = A_laplace(datos,C); %Pseudo-pressure matrix A. [Vx*Vy, Vx*Vy]
matriu_A(1,1)=-5;

pseudo_p = inv(matriu_A)*div_u_p; %pseudo_p = [Vx*Vy,1]

[p_gradX, p_gradY] = gradient_p(datos, pseudo_p, nodal_mesh); 
p_gradX = haloupdate(p_gradX);
p_gradY = haloupdate(p_gradY);
u_p = haloupdate(u_p);
v_p = haloupdate(v_p);

u_n1 = u_p - p_gradX;
v_n1 = v_p - p_gradY;

divergence_field = divergencia_u(datos, u_n1, v_n1, nodal_mesh, num);
% 
sum(divergence_field)
    
   