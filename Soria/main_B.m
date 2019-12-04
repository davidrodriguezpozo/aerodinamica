function main_B

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


Vx = 3; %nombre de divisions en x de V.C SENSE HALO
Vy = 3;



datos = INPUT(Vx,Vy);
C = meshes (datos, datos.malla);

matriu_A = A_laplace(datos,C); %Pseudo-pressure matrix A.

v = zeros(Vx,Vy);
u = zeros(Vx,Vy);
u(2,2) = 1;

[nodal_mesh num] = nodalmesh(Vx,Vy)

u_p = divergencia_u(datos, u, v, nodal_mesh);

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

for i = 1:datos.Nx
    for j = 1:datos.Ny
        
        x = C.stagX_x(i,j);
        y = C.stagX_y(i,j);
        
        u(i,j) = datos.F*cos(2*pi*x)*sin(2*pi*y);
        
        x = C.stagY_x(i,j);
        y = C.stagY_y(i,j);
        
        v(i,j) = -datos.F*cos(2*pi*y)*sin(2*pi*x);
        
    end
end

function [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v)

        conv_u = convection_main(u,v,datos.L,datos.H);
        diff_u = diffusion_u(u,datos.L);

        conv_v = convection_main(v,u,datos.H,datos.L);
        diff_v = diffusion_u(v,datos.H);


function cu = convection_main (u,v,L,H)
    
        Nx = size(u,1);
        Ny = size(u,2);
        
        dx = L/Nx; dy = H/Ny;

    for i = 2:Nx-1
        for j=2:Ny-1

        % CDS : Central Difference Scheme
        ue = (u(i+1,j)+u(i,j))/2;
        uw = (u(i+1,j)+u(i,j))/2;
        vn = (v(i,j+1)+u(i,j))/2;
        vs = (v(i,j-1)+u(i,j))/2;
        un = (u(i,j+1)+u(i,j))/2;
        us = (u(i,j-1)+u(i,j))/2;
   
        % Mass Fluxes 
   
        Fe = ue * dy;
        Fw = uw * dy;
        Fn = vn * dx;
        Fs = vs * dx;
    
        % Convection term
   
        cu(i-1,j-1) = Fe*ue - Fw*uw + Fn*un - Fs*us;
        end
    end

    
function [conv diff] = Analytic (datos, C, malla)
% This function ... wh
    for i = 2:datos.Nx-1
        for j=2:datos.Ny-1
        
        dx = C.dx(i);
        dy = C.dy(j);
        Sup = dx*dy;
        
        if strcmp(malla,'x') == 1
            x = C.stagX_x(i,j);
            y = C.stagX_y(i,j);
        elseif strcmp(malla,'y') == 1
            x = C.stagY_x(i,j);
            y = C.stagY_y(i,j);
        end
        
        u = datos.F*cos(2*pi*x)*sin(2*pi*y);
        v = -datos.F*cos(2*pi*y)*sin(2*pi*x);
        
        du_dx = - datos.F*2*pi*sin(2*pi*x)*sin(2*pi*y);
        du_dy = datos.F*2*pi*cos(2*pi*x)*cos(2*pi*y);
        dv_dx = -datos.F*2*pi*cos(2*pi*y)*cos(2*pi*x);
        dv_dy = datos.F*2*pi*sin(2*pi*y)*sin(2*pi*x);
        
        d2u_dx2 = - datos.F^2*4*pi^2*sin(2*pi*y)*cos(2*pi*x);
        d2v_dy2 = datos.F^2*4*pi^2*sin(2*pi*x)*cos(2*pi*y);
        
        if strcmp(malla,'x') == 1
            
            conv(i-1,j-1) = Sup*(u*du_dx + u*dv_dx); %u * du/dx 
            diff(i-1,j-1) = Sup*(d2u_dx2); 
            
        elseif strcmp(malla,'y') == 1
            
            conv(i-1,j-1) = Sup*(v*du_dy + v*dv_dy); %v * du/dx 
            diff(i-1,j-1) = Sup*(d2v_dy2); 
        
        end 
    end
    end


function datos = INPUT(Vx,Vy)

    datos.Nx = Vx+2;
    datos.Ny = Vy+2;
    
    datos.Vx = Vx;
    datos.Vy = Vy;
    
    datos.L = 1;
    datos.H = 1;
    
    datos.uniform = true;
    datos.gamma = 1;
    
    datos.malla = 2;
    
    datos.F = 1;
    
    
function error = ERROR(numeric,analytic)
    
    error = 0;

    for i = 1:size(numeric,1)
       for j = 1:size(numeric,2)
           %if analytic(i,j) < 0.0001 ||  analytic(i,j) > -0.0001
            %    a = abs(analytic(i,j)-numeric(i,j));
             %   disp('Algo va mal');
           %else 
               a = abs(analytic(i,j)-numeric(i,j));
           %end
            if a > error 
                error = a;
            end
       end
    end
   
   