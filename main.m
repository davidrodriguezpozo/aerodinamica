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

datos = INPUT(Vx,Vx);
C = meshes (datos, datos.malla);

matriu_A = zeros(datos.Nx*datos.Nx,datos.Nx*datos.Nx);
matriu_A = A_laplace(datos,C);


% FEM SOLUCIO ANALITICA
malla = 'x';
[S.cu_anal S.du_anal] = Analytic (datos, C, malla);

malla = 'y';
[S.cv_anal S.dv_anal] = Analytic (datos, C, malla);

%disp('Analitica calculada');

clear malla

% PERFIL VELOCITATS

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

%A = Poisson (datos, C)

if k == length(div)
    figure
    contourf(C.stagY_x,C.stagY_y,v,'LineWidth',0.1) ;
    c = colorbar;
    str = {'Y-Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    
    figure
    contourf(C.stagX_x,C.stagX_y,u,'LineWidth',0.1) ;
    c = colorbar;
    str = {'X-Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1
            x_v (i-1,j-1) = C.X(i);
            y_v (i-1,j-1) = C.Y(j);
            
            vel_x = (u(i,j) + u(i-1,j))/2;
            vel_y = (v(i,j) + v(i,j-1))/2;
            
            modVel (i-1,j-1) = sqrt(vel_x^2 + vel_y^2);
            
        end
    end
    
    figure
    contourf(x_v,y_v,modVel,'LineWidth',0.1) ;
    c = colorbar;
    str = {'Module of Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
end

% FEM SOLUCIO NUMERICAL

[S.cu S.diffu S.cv S.diffv] = Numerical (datos, C, u, v);

%disp('Numerica calculada');


error_mat_cu = abs(S.cu - S.cu_anal);

error_du(k) = ERROR(S.diffu,S.du_anal);
error_cu(k) = ERROR(S.cu,S.cu_anal);
error_dv(k) = ERROR(S.diffv,S.dv_anal);
error_cv(k) = ERROR(S.cv,S.cv_anal);
%Nx_vector(i) = datos.Nx;
%disp('Error calculat');
end

figure;
loglog(1./Nx_vector,error_du,'o-');
title('Diffusive term error (x-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error diffusive term (x-dir)');

figure;
loglog(1./Nx_vector,error_dv,'o-');
title('Diffusive term error (y-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error diffusive term (y-dir)');

figure;
loglog(1./Nx_vector,error_cu,'o-');
title('Convective term error (x-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error convective term (x-dir)');

figure;
loglog(1./Nx_vector,error_cv,'o-');
title('Convective term error (y-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error convective term (y-dir)');





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
   
   