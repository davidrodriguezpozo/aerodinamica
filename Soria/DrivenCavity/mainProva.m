%% DRIVEN CAVITY PROBLEM 

% Code developed by:
% - Sergi Martínez Castellarnau
% - Carlos Pérez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

% Execute function main 

function main

clear all
clc
close all


datos = INPUT;

C = MESH (datos);

S = CondicionesIniciales(datos,C);

[ae aw an as ap bp] = Coefficients (datos);

S = Solve_NS(datos,C, S);


function [R_u R_v R_u_prev R_v_prev] = factor_R  (datos, S)
    
    dx =datos.L/datos.Nx;
    
    % Update 
    
    R_u_prev = S.R_u; 
    R_v_prev = S.R_v; 
    
    for i = 2:datos.Nrow_sx -1
        for j = 2:datos.Ncol_sx -1
            %R_u(i,j) = conv(u) - diff(u);
            R_u(i,j) = 1;
        end
    end
    
    for i = 2:datos.Nrow_sy -1
        for j = 2:datos.Ncol_sy -1
            %R_v(i,j) = conv(u) - diff(u);
            R_v(i,j) = 1;
        end
    end
  
    %{
    for (i=1; i<=N_row_stag_x-2; i++){
        for (j=1; j<=N_col_stag_x-2; j++){
            R_U[i][j]=mu*(U[i][j+1]+U[i][j-1]+U[i+1][j]+U[i-1][j]-4*U[i][j])- rho*0.25*dx*((U[i][j]+U[i][j+1])*(U[i][j]+U[i][j+1])-(U[i][j]+U[i][j-1])*(U[i][j]+U[i][j-1])+(U[i][j]+U[i-1][j])*(V[i][j-1]+V[i][j])- (U[i][j]+U[i+1][j])*(V[i+1][j-1]+V[i+1][j]));
        }
    }

    for (i=1;i<=N_row_stag_y-2;i++){
        for (j=1;j<=N_col_stag_y-2;j++){
            R_V[i][j]=mu*(V[i][j+1]+V[i][j-1]+V[i-1][j]+V[i+1][j]-4*V[i][j])-rho*0.25*dx*((U[i][j+1]+U[i-1][j+1])*(V[i][j]+V[i][j+1])-(U[i][j]+U[i-1][j])*(V[i][j]+V[i][j-1])+(V[i][j]+V[i-1][j])*(V[i][j]+V[i-1][j])- (V[i][j]+V[i+1][j])*(V[i][j]+V[i+1][j]));
        }
    }
    %}


function datos = INPUT

datos.Re = 100;

datos.Nx = 5;
datos.Ny = 5;

datos.Ncol = datos.Nx + 2; % Numero de nodos en x-direction
datos.Nrow = datos.Ny + 2; % Numero de nodos en x-direction

datos.Ncol_sx = datos.Ncol - 1;
datos.Nrow_sx = datos.Nrow - 2;

datos.Ncol_sy = datos.Ncol - 2;
datos.Nrow_sy = datos.Nrow - 1;

datos.L = 1;
datos.H = 1;

datos.tfinal = 100;

datos.dt_inicial = 0.01;
datos.rho = 1;
datos.mu = 1/datos.Re;


function C = MESH (datos)

fprintf('-( 1 )- Calculando mallas...\n');

% Preallocation of variables

C.X = zeros (datos.Nrow,datos.Ncol);
C.Y = zeros (datos.Nrow,datos.Ncol);

C.dx = zeros (datos.Nrow,1);
C.dy = zeros (datos.Ncol,1);

C.stagX_x = zeros (datos.Nrow_sx,datos.Ncol_sx);
C.stagX_y = zeros (datos.Nrow_sx,datos.Ncol_sx);

C.stagY_x = zeros (datos.Nrow_sy,datos.Ncol_sy);
C.stagY_y = zeros (datos.Nrow_sy,datos.Ncol_sy);

[ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( datos );

[ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( datos, C );

[ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( datos, C );

fprintf('---> Mallas calculada\n')


function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( datos )

% Delta X

for i = 1:datos.Ncol

dx(i) = spacing_no_uniform( datos.L, i, datos.Nx ) - spacing_no_uniform( datos.L, i-1, datos.Nx );

end

% Delta Y

for j = 1:datos.Nrow

dy(j) = spacing_no_uniform( datos.H, j, datos.Ny ) - spacing_no_uniform(datos.H, j-1, datos.Ny );

end

% COLLOCATED MESHES

X (1) = 0; X(2) = dx(1)/2;

for i = 3:datos.Ncol-1
    X(i) = X(i-1) + dx(i-1);
end
X(datos.Ncol) = X(datos.Ncol-1) + dx(datos.Ncol)/2;

Y (1) = 0; Y(2) = dy(1)/2;
for j = 3:datos.Nrow-1
    Y(j) = Y(j-1) + dy(j-1);
end
Y(datos.Nrow) = Y(datos.Nrow-1) + dy(datos.Nrow)/2;

for i=1:datos.Nrow
    for j=1:datos.Ncol
        Coll_X(i,j) = X(i);
        Coll_Y(i,j) = Y(j);
    end
end


function [ X, Y ] = STAGG_MESH_X ( datos, C )

vector_x = zeros( datos.Nx+1 ,1 );
vector_x(1) = 0; % La malla comença a la paret
for i = 2:datos.Nx
    vector_x(i) = vector_x(i-1) + C.dx(i);
end
vector_x(datos.Nx+1) = datos.L;

for i = 1:datos.Nrow_sx+1
    for j = 1:datos.Ncol_sx
        X(i,j) = vector_x(j);
        Y(i,j) = C.Y(i);
    end
end

fprintf(' Staggered Mesh X calculada\n');

function [ X, Y ] = STAGG_MESH_Y ( datos, C )

vector_y = zeros( datos.Ncol_sy ,1 );
vector_y(1) = 0; % La malla comença a la paret

for i = 2:datos.Ny
    vector_y(i) = vector_y(i-1) + C.dx(i);
end
vector_y(datos.Ny+1) = datos.L;

for i = 1:datos.Nrow_sy
    for j = 1:datos.Ncol_sy+1
        X(i,j) = C.X(j);
        Y(i,j) = vector_y(i);
    end
end

fprintf(' Staggered Mesh Y calculada\n');

function delta = spacing_no_uniform (L, i, N )
    delta = L * i / N;


function S = CondicionesIniciales (datos, C)

fprintf('\n-( 2 )- Condiciones iniciales establecidas \n');

S.P = zeros(datos.Nrow, datos.Ncol);
S.U = zeros(datos.Nrow_sx, datos.Ncol_sx);
S.V = zeros(datos.Nrow_sy, datos.Ncol_sy);
S.P_prev = ones(datos.Nrow, datos.Ncol);
S.U_prev = zeros(datos.Nrow_sx, datos.Ncol_sx);
S.V_prev = zeros(datos.Nrow_sy, datos.Ncol_sy);
S.R_u = zeros(datos.Nrow_sx, datos.Ncol_sx);
S.R_v = zeros(datos.Nrow_sy, datos.Ncol_sy);
S.R_u_prev = zeros(datos.Nrow_sx, datos.Ncol_sx);
S.R_v_prev = zeros(datos.Nrow_sy, datos.Ncol_sy);
S.time = 0;
S.dt = datos.dt_inicial;


function S = Solve_NS (datos, C, S)

fprintf('\n-( 3 )- Calculando solución... \n');

time = datos.dt_inicial;

while max_dif < datos.delta || time < datos.tfinal 
    
    [S.R_U S.R_V S.R_U_prev S.R_V_prev] = factor_R  (datos, S);

    time = CalculateTimeStep (datos, S.U, S.V, time);
    time = 1000;
end

disp('----> Solución calculada \n')


function [ae aw an as ap bp] = Coefficients (datos)
    
    for i=1:datos.Nrow
        for j=1:datos.Ncol
            ae(i,j)=1; aw(i,j)=1;
            an(i,j)=1; as(i,j)=1;
            bp(i,j)=0; 
        end
    end
    
    % BOUNDARIES
    
    for i = 1:datos.Nrow
        aw(i,1) = 0;
        an(i,1) = 0;
        as(i,1) = 0;

        ae(i,datos.Ncol) = 0;
        an(i,datos.Ncol) = 0;
        as(i,datos.Ncol) = 0;
    end
    
    for i=2:datos.Nrow-2 
        ae(i,1)=2; aw(i,1)=2;
        aw(i,datos.Ncol)=2; ae(i,datos.Ncol-1)=2;
    end
        
    for j = 1:datos.Ncol
        an(1,j) = 0;
        an(1,j) = 0;
        as(1,j) = 0;

        ae(datos.Nrow,j) = 0;
        an(datos.Nrow,j) = 0;
        as(datos.Nrow,j) = 0;
    end
    
    for j = 2:datos.Ncol-1
        as(1,j)=2; an(1,j)=2;
        an(datos.Nrow,j)=2; as(datos.Nrow-1,j)=2;
    end
    
    for i=1:datos.Nrow
        for j = 1: datos.Ncol
            ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j);
        end
    end
    
    


%{

while converg == false || S.time < datos.tfinal

% Update de variables

S.U_prev = S.U;

S.V_prev = S.V;

% Factor R

[S.R_u s.R_v] = Compute_R(datos,C,S);

end

%}

% Print results

function [R_u R_v] = Compute_R (datos,C,S)

% Update de variables

S.R_u_prev = S.R_u;

S.R_v_prev = S.R_v;

%% NOMENCLATURA

% Fe = flux sortint massic en la cara e (EAST)

% uE = velocitat direccio x en el node E (EAST)

% ue = velocitat direccio x en la cara e (EAST)

% Se = superficie de la cara e (EAST)

% Vegeu que sub index e indica cara; E indica node

for i=1:datos.Nrow_sx

for j=1:datos.Ncol_sx

UE = S.U(i,j+1);

UW = S.U(i,j-1);

UN = S.U(i+1,j);

US = S.U(i-1,j);

UP = S.U(i,j+1);

Se = C.stagX_y(i,j+1) - C.stagX_y(i,j);

Sw = C.stagX_y(i,j) - C.stagX_y(i,j);

Sn = C.stagX_x(i,j) - C.stagX_x(i,j);

Ss = C.stagX_x(i,j) - C.stagX_x(i,j);

% Aproximacio terme convectiu

ue = (U(i,j+1)+U(i,j))/2; Fe = ue*Se;

uw = (U(i,j+1)+U(i,j))/2; Fw = uw*Sw

un = (U(i,j+1)+U(i,j))/2;

us = (U(i,j+1)+U(i,j))/2;

vn = (V(i,j+1)+V(i,j))/2; Fn = vn*Sn;

vs = (V(i,j+1)+V(i,j))/2; Fs = vs*Ss;

end

end


function time = CalculateTimeStep (datos, u, v, time)

    dx = datos.L/datos.Nx;
    
    Conv_dt = 0; Diff_dt = 0;
    
    umax = 0; vmax = 0;
    
    for i = 1:datos.Nrow_sx
        for j = 1:datos.Ncol_sx
            if u(i,j) > umax
                umax = u(i,j);
            end
        end
    end
    
    for i = 1:datos.Nrow_sy
        for j = 1:datos.Ncol_sy
            if v(i,j) > vmax
                vmax = v(i,j);
            end
        end
    end
    
    modVel = sqrt (umax*umax + vmax*vmax);
    
    Conv_dt = 0.35*dx/modVel;
    Diff_dt = 0.2*datos.rho*dx*dx/datos.mu; 
    
    dt = min(Conv_dt ,Diff_dt)*0.8; 
    
    time = time + dt;
