function [p pseudo_p delta_t] = pressure (datos, u_p, v_p)

% Vx = 3; %nombre de divisions en x de V.C SENSE HALO
% Vy = 3;
% R_uant = zeros(Vx+2, Vx+2);
% R_vant = zeros(Vx+2, Vx+2);
% R_u = zeros(Vx+2, Vx+2);
% R_v = zeros(Vx+2, Vx+2);
% Re = 10000;
% datos = INPUT(Vx,Vy,Re);
% C = meshes (datos, datos.malla);
% 
% v = zeros(datos.Nx,datos.Ny);
% u = zeros(datos.Nx,datos.Ny);
% u(3,3) = 1; % El camp de velocitats t� HALO
% v(3,3) = 1;
% u = haloupdate(u);
[nodal_mesh num] = nodalmesh(datos.Vx,datos.Vy);
%[nodal_mesh num] = nodalmesh(datos.Vx,datos.Vy);
delta_t = TimeStep(datos,u_p,v_p);

div_u_p = divergencia_u(datos, u_p, v_p, nodal_mesh, num); %u_p = [Vx*Vy,1]

pseudo_p = zeros(datos.Vx*datos.Vx,1); %Pseudo-pressure vector NOT KNOWN 


matriu_A = A_laplace(datos); %Pseudo-pressure matrix A. [Vx*Vy, Vx*Vy]
matriu_A(1,1)=-5; %Per solucionar el problema de la matriu A singular (no es podria invertir)

pseudo_p = inv(matriu_A)*div_u_p; %pseudo_p = [Vx*Vy,1]

p = pseudo_p*datos.rho / delta_t;

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
    
   