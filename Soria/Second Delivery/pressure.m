function [P_matrix_halo, pseudo_p, delta_t] = pressure (datos, u_p, v_p, nodal_mesh, delta_t)

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
u(3,3) = 1; % El camp de velocitats t� HALO
% v(3,3) = 1;
% u = haloupdate(u);
[nodal_mesh, num] = nodalmesh(datos.Vx,datos.Vy);
%[nodal_mesh num] = nodalmesh(datos.Vx,datos.Vy);
delta_t = TimeStep(datos,u_p,v_p,delta_t);

div_u_p = divergencia_u(datos, u_p, v_p, nodal_mesh, num); %u_p = [Vx*Vy,1]

pseudo_p = zeros(datos.Vx*datos.Vx,1); %Pseudo-pressure vector NOT KNOWN 


matriu_A = A_laplace(datos); %Pseudo-pressure matrix A. [Vx*Vy, Vx*Vy]
matriu_A(1,1)=-5; %Per solucionar el problema de la matriu A singular (no es podria invertir)

pseudo_p = matriu_A\div_u_p; %pseudo_p = [Vx*Vy,1]

p = pseudo_p*datos.rho / delta_t;

P_matrix = zeros(datos.Vx,datos.Vy);

for i=1:datos.Vy
    for j=1:datos.Vx
        k = nodal_mesh(i,j);
        P_matrix(i,j) = p(k); % P_matrix =  [Vx,Vy]
    end
end

P_matrix_halo = zeros(datos.Vx+2,datos.Vy+2);

for i=2:datos.Vx+1 
    for j=2:datos.Vy+1
    
    P_matrix_halo(i,j) = P_matrix(i-1,j-1);

    end
end

    
   