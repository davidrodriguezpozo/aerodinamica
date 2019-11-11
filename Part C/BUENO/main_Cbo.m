function main_C

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');
set(groot, 'DefaultLegendInterpreter','latex');

clear all 
clc
close all 


Vx = 10; 
Vy = 10;
Re = 1000;

datos = INPUT(Vx,Vy,Re);
C = meshes (datos, datos.malla);

matriu_A = A_laplace(datos,C); %Pseudo-pressure matrix A.

time = 0.00;
datos.F = exp(-8*pi^2*datos.mu*time);

u = zeros(datos.Nx, datos.Ny);
v = zeros(datos.Nx, datos.Ny);
for i = 2:datos.Nx-1
    for j = 2:datos.Ny-1
        
        x1 = C.stagX_x(i,j);
        y1 = C.stagX_y(i,j);
        
        u(i,j) = datos.F*cos(2*pi*x1)*sin(2*pi*y1);
        
        x2 = C.stagY_x(i,j);
        y2 = C.stagY_y(i,j);
        
        v(i,j) = -datos.F*cos(2*pi*y2)*sin(2*pi*x2);
        
    end
end

u = haloupdate(u);
v = haloupdate(v);


%STEP TIME FOR STABILITY:
%nu = 1.5e-5 -> Viscosidad cinemática [m2/s] aire
nu=0.010;
datos.mu = 0.010;

time = 0;

final_time = 2;

R_u = zeros(datos.Nx);
R_v = zeros(datos.Nx);
p = zeros(Vx^2,1);
p_prev = zeros(Vx^2,1);


conv_u = zeros(datos.Nx); 
diff_u = zeros(datos.Nx);
conv_v = zeros(datos.Nx);
diff_v = zeros(datos.Nx);

u_p = zeros(datos.Nx);
v_p = zeros(datos.Nx);
u_prev = zeros(datos.Nx);
v_prev = zeros(datos.Nx);

ae = zeros(datos.Nx);
aw = zeros(datos.Nx);
an = zeros(datos.Nx);
as = zeros(datos.Nx);
ap = zeros(datos.Nx);
bp = zeros(datos.Nx);

P = zeros(datos.Nx);
P_prev = zeros(datos.Nx);
dif_P = 1;
dif_U = 1;
dif_V = 1;
delta_V = 1e-7;

Nx = datos.Nx;
Ny = datos.Ny;

Td = 0.5*(datos.L/datos.Vx)*(datos.H/datos.Vy) / datos.mu;
delta_T = 0.02;
time = delta_T;
    
index = 1;

while time <= final_time && dif_U > delta_V && dif_V > delta_V
    
    R_u_prev = R_u;
    R_v_prev = R_v;
    
    u_prev = u;
    v_prev = v;
    
    %datos.F = exp(-4*pi*datos.mu/datos.rho*time);
    
    [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);
    
    %Calculem R: R = - Conv(u)/V + Diff(u)/V
    Vol = (datos.L/datos.Vx)*(datos.H/datos.Vy);
    R_u = -(conv_u/Vol) + datos.mu*(diff_u/Vol);
    R_v = -(conv_v/Vol) + datos.mu*(diff_v/Vol);
    
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1
        u_p(i,j) = u(i,j) + delta_T*(3/2*R_u(i,j) - 1/2*R_u_prev(i,j));
        v_p(i,j) = v(i,j) + delta_T*(3/2*R_v(i,j) - 1/2*R_v_prev(i,j));
        end
    end
    
    u_p = haloupdate(u_p);
    v_p = haloupdate(v_p);
    
    dx = datos.L/datos.Vx;
    dy = datos.H/datos.Vy;
    
    %delta_T = 0.02;
    
    rho = 1;
    
    for i=2:Nx-1
        for j=2:Ny-1
            ae(i,j) = dy/dx;
            aw(i,j) = dy/dx;
            an(i,j) = dx/dy;
            as(i,j) = dx/dy; 
            bp(i,j)=-1/delta_T*rho*(u_p(i,j)*(dy)-u_p(i,j-1)*(dy)+v_p(i,j)*(dx)-v_p(i-1,j)*(dx));
            
        end
    end
    
    ap = ae + aw + an + as;
   
    
    [P P_prev dif_P] = SolverP (ae, aw, an, as, ap, bp, P_prev, P, dif_P, datos);
    
    P = haloupdate(P);
    
    p_analytic = Analytic_Pressure (datos, C, time);
    
    rho = datos.rho;
    mu = datos.mu;
   
    index = index+1;
    t_vector(index)=time;
    u_vector(index)=u(3,3);
    v_vector(index)=v(3,3);
    P_anal(index)=-(exp(-8*pi^2*mu*time))^2*(rho*(cos(4*pi*C.Coll_x(3,3))+cos(4*pi*C.Coll_y(3,3))))/4;
    P_vec(index)=P(3,3);
    u_anal(index) = datos.F*cos(2*pi*C.stagX_x(3,3))*sin(2*pi*C.stagX_y(3,3));
    
    errorP(index) = abs(P_anal(index) - P_vec(index));
    acu_time (index) = time;
    
    for i=2:Nx-1
        for j=2:Nx-1
 
        u(i,j)=u_p(i,j)-delta_T/rho*(P(i,j+1)-P(i,j))/(dx);
        v(i,j)=v_p(i,j)-delta_T/rho*(P(i+1,j)-P(i,j))/(dy);
 
        end
    end
    
    u = haloupdate(u);
    v = haloupdate(v);
    
    dif_U = 0;
    dif_V = 0;
    
    for i = 1:Nx
        for j = 1:Ny
            if abs(u(i,j)-u_prev(i,j))> dif_U
                dif_U = abs(u(i,j)-u_prev(i,j));
            end
            
            if abs(v(i,j)-v_prev(i,j))> dif_V
                dif_V = abs(v(i,j)-v_prev(i,j));
            end
        end
    end
    
        %Time for convective term
    Tc_u = min(datos.L/(datos.Vx*max(max(u))));
    Tc_v = min(datos.H/(datos.Vy*max(max(v))));
    Tc = min (Tc_u,Tc_v);
    %Time for diffusive term
    Td = 0.5*(datos.L/datos.Vx)*(datos.H/datos.Vy) / datos.mu;
    
    %Step time:
    delta_T =0.2 * min(Tc,Td);
    delta_T =0.02 ;
    
    %delta_T = 0.0421; %per 5
    
    time = time + delta_T;
    
    datos.F = exp(-8*pi^2*datos.mu*time);
    

end

% 
semilogy(acu_time,P_anal); hold on;
semilogy(acu_time,P_vec); 
legend('Analytic','Numerical')

figure; 
semilogy(acu_time,u_anal); hold on;
semilogy(acu_time,u_vector);
legend('Analytic','Numerical');
posx = num2str(C.stagX_x(3,3));
posy = num2str(C.stagX_y(3,3));
title(['Horizontal Velocity at point x =',posx,' and y = ',posy]);

disp('AA');
figure;
    contourf(C.Coll_x,C.Coll_y,P,'LineWidth',0.05) ;
    c = colorbar;
    str = {'Pressure'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Pressure','Fontsize',16)



function A = A_laplace(datos,C)

Vx = datos.Nx;
Vy = datos.Ny;

[nodal_mesh num] = nodalmesh(Vx,Vy);


   
    A = zeros(num,num); %Serà la matriu que es multiplicarà per p.
    
    for i=1:num
        %Diagonal term:
        A(i,i) = -4;
        
       [k j] = find(nodal_mesh == i);     
       %NORD:
       if k-1 == 0 %Aquest node es troba a la frontera esquerra
           A(i, nodal_mesh(Vy,j)) = 1;
       else
           A(i, nodal_mesh(k-1,j)) = 1;
       end

       %SUD:
       if k+1 > Vy
           A(i, nodal_mesh(1,j)) = 1;        
       else
           A(i, nodal_mesh(k+1,j)) = 1;
       end

       %WEST:
       if j-1 == 0
           A(i, nodal_mesh(k,Vx)) = 1;
       else
           A(i, nodal_mesh(k,j-1)) = 1;
       end

       %EAST:
       if j+1 > Vx            
           A(i, nodal_mesh(k,1)) = 1;
       else
           A(i, nodal_mesh(k,j+1)) = 1;
       end
    end



    
function [nodal_mesh num] = nodalmesh(Vx,Vy)
%{
this function creates the matrix "nodal_mesh", which contains the indexes
%of each node, ordered as follows: 
    
              7  8  9  
nodal_mesh = [4  5  6]
              1  2  3
%}

    nodal_mesh = zeros(Vx,Vy); %Matriu de tots els nodes dels VC.
    num = Vx*Vy; %Per sabre fins a quin index arribem
    ind=1;
    for i=Vx:-1:1
        for j=1:Vy
            nodal_mesh(i,j)=ind; %S'ompla la matriu de nodes amb els index necessaris;
            ind = ind+1;
        end
    end



