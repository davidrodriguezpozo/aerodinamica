function main

clear all 
clc
close all 

datos = INPUT;
C = meshes (datos, datos.malla);

% FEM SOLUCIO ANALITICA
malla = 'x';
[S.cu_anal S.du_anal] = Analytic (datos, C, malla);

malla = 'y';
[S.cv_anal S.dv_anal] = Analytic (datos, C, malla);

disp('Analitica calculada');

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

% FEM SOLUCIO NUMERICAL

[S.cu S.diffu S.cv S.diffv] = Numerical (datos, C, u, v);

disp('Numerica calculada');

error_du = ERROR(S.diffu,S.du_anal);
error_cu = ERROR(S.cu,S.cu_anal);
error_dv = ERROR(S.diffv,S.dv_anal);
error_cv = ERROR(S.cv,S.cv_anal);

disp('Error calculat');

function [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v)

        conv_u = convection(u,v,datos.L,datos.H);
        diff_u = diffusion_u(u,datos.L);

        conv_v = convection(v,u,datos.H,datos.L);
        diff_v = diffusion_u(v,datos.H);



function cu = convection (u,v,L,H)
    
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

    for i = 2:datos.Nx-1
        for j=2:datos.Ny-1
        
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
        du_dy = datos.F*2*pi*sin(2*pi*x)*cos(2*pi*y);
        dv_dx = -datos.F*2*pi*sin(2*pi*y)*sin(2*pi*x);
        dv_dy = datos.F*2*pi*sin(2*pi*y)*sin(2*pi*x);
        
        d2u_dx2 = - datos.F^2*4*pi^2*sin(2*pi*y)*cos(2*pi*x);
        d2v_dy2 = datos.F^2*4*pi^2*sin(2*pi*x)*cos(2*pi*y);
        
        if strcmp(malla,'x') == 1
            
            conv(i-1,j-1) = u*du_dx + u*dv_dx; %u * du/dx 
            diff(i-1,j-1) = d2u_dx2; 
            
        elseif strcmp(malla,'y') == 1
            
            conv(i-1,j-1) = v*du_dy + v*dv_dy; %v * du/dx 
            diff(i-1,j-1) = d2v_dy2; 
        
        end 
    end
end

function datos = INPUT 
    datos.Vx = 3;
    datos.Vy = 3;
    datos.Nx = datos.Vx+2;
    datos.Ny = datos.Vy+2;
    
    datos.L = 1;
    datos.H = 1;
    
    datos.uniform = true;
    datos.gamma = 1;
    
    datos.malla = 3;
    
    datos.F = 1;
    
    
function error = ERROR(numeric,analytic)
    
    error = 0;

    for i = 1:size(numeric,1)
       for j = 1:size(numeric,2)
            a = abs(analytic(i,j)-numeric(i,j))/analytic(i,j); 
            
            if a > error 
                error = a;
            end
       end
    end
   