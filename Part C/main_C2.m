function main_C

%Es la mateixa funcio que el main pero sense les divisions de VC pel M.M.S.

%{  
                                El Codi hauria de fer aix? 

 1. Timestep for stability
 2. u_n_1 i u_n conegudes --> R_n_1 i R_n
 3. Evaluar u_p = u_n + ?t(3/2 * R_n - 1/2 * R_N_1
        R = - Conv(u)/V + Diff(u)/V
 4. Resoldre l'equaci? de Poisson A * pseudo_p = Div(u_p)
 5. Obtenir pseudo_p
 6. Obtenir u_n1 --> u_n1 = u_p - grad(pseudo_p) 
 7. timestep + 1 
 %}

clear all 
clc
close all 


Vx = 20; %N? de divisions en x de V.C.
Vy = 20;

datos = INPUT(Vx,Vy);
C = meshes (datos, datos.malla);

matriu_A = A_laplace(datos,C); %Pseudo-pressure matrix A.

%Perfil de velocitats:
% PERFIL VELOCITATS

time = 0;
datos.F = exp(-8*pi^2*datos.mu*time);

u = zeros(datos.Nx, datos.Ny);
v = zeros(datos.Nx, datos.Ny);
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


%STEP TIME FOR STABILITY:
%nu = 1.5e-5 -> Viscosidad cinemática [m2/s] aire
nu=1.5e-5;

time = 0;

final_time = 10;

R_u = zeros(Vx);
R_v = zeros(Vx);
p = zeros(Vx^2,1);
p_prev = zeros(Vx^2,1);

while time <= final_time 
    
    R_u_prev = R_u;
    R_v_prev = R_v;
    
    %Time for convective term
    %Tc = min((datos.L/(datos.Nx*max(max(u)))), (datos.H/(datos.Vy*max(max(v)))));
    Tc = min((datos.L/(datos.Nx*max(max(u)))), (datos.H/(datos.Vy*max(max(v)))));
    %Time for diffusive term
    Td = 0.5*(datos.L/datos.Nx)*(datos.H/datos.Ny) / nu;

    %Step time:
    delta_T =0.4 * min(Tc,Td);
    
    time = time + delta_T;
    
    datos.F = exp(-8*pi^2*datos.mu*time);
    
    %Calculem R: R = - Conv(u)/V + Diff(u)/V
    [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);

    Vol = (datos.L/datos.Vx)*(datos.H/datos.Vy);
    R_u = -(conv_u/Vol) + datos.mu*(diff_u/Vol);
    R_v = -(conv_v/Vol) + datos.mu*(diff_v/Vol);

    [nodal_mesh num] = nodalmesh(Vx,Vy);
    
    u_p = u(2:end-1,2:end-1) + delta_T*(3/2*R_u - 1/2*R_u_prev);
    v_p = v(2:end-1,2:end-1) + delta_T*(3/2*R_v - 1/2*R_v_prev);
    
    [nodal_mesh num] = nodalmesh(Vx,Vy);
    
    div_u_p = divergencia_u(datos, u_p, v_p, nodal_mesh, num);
    
    sum(div_u_p)
    
    pseudo_p = zeros(Vx*Vx,1); %Pseudo-pressure vector NOT KNOWN 

    %Per solucionar el problema de la matriu A singular (no es podria invertir)

    matriu_A(1,1)=-5;
    
    %pseudo_p = inv(matriu_A)*u_p;
    pseudo_p = inv(matriu_A)*div_u_p;
    
    rho = 1;
    
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
    
    p_analytic = Analytic_Pressure (datos, C);
    
    error = ERROR(p_mat,p_analytic)

    [p_gradX p_gradY] = gradient_p(datos, pseudo_p, nodal_mesh, num);
    
    u_n1 = u(2:end-1,2:end-1) - p_gradX;
    v_n1 = v(2:end-1,2:end-1) - p_gradY;
    
    %% UPDATE OF VELOCITY FIELD
    
    u(2:end-1,2:end-1) = u_n1;
    v(2:end-1,2:end-1) = v_n1;
    
    divvvvvv = divergencia_u(datos, u_n1, v_n1, nodal_mesh, num);

    sum(divvvvvv)

end

disp('AA');



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
 
function A = A_laplace(datos,C)

Vx = datos.Vx;
Vy = datos.Vy;

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


function p = Analytic_Pressure (datos, C)
% This function ... wh
    for i = 2:datos.Nx-1
        for j=2:datos.Ny-1
        
        x = C.X(i);    
        y = C.Y(j);
        p(i-1,j-1) = - datos.F^2 * datos.rho * (cos(4*pi*x)+ cos(4*pi*y))/4;
           
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
    
    datos.mu = 1.5e-5;
    datos.rho = 1;
    
    
    
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION MESHES         %%%%%%%%%%%%%
%%% This function allows the programer to allocate the  %%%
%%% different points of the three meshes needed to sol- %%%
%%% ve the NS equations. Furthermore, three functions   %%%
%%% are contempled in order to fullfill all the require-%%%
%%% ments. First of all, a Colloated-Mesh function will %%% 
%%% be performed moving the y-axis center and the x-axis%%%
%%% center up and right, respectively. Afterwards, two  %%%
%%% Staggered-Mesh functions will be executed (from x & %%%
%%% y-axis) building two other meshes into a "regular"  %%%
%%% performing position (they are explained below). Fin-%%%
%%% ally a position-incremental function is defined and %%%
%%% the sum of all is explained in the meshes function, %%%
%%% that works as a "main".                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  C = meshes (datos, malla)
    % A C class is defined in order to organise all the outputs of the
    % function.
  
        
    % Preallocation of variables: X and Y positions (corresponding to
    % collocated locations), incrementals, staggered (x and y) matrixes
    % which will allocate the nodal positions.
    C.X = zeros (datos.Nx,datos.Ny);
    C.Y = zeros (datos.Nx,datos.Ny);    
    C.dx = zeros (datos.Nx,1);
    C.dy = zeros (datos.Ny,1);
    
    % Inicializtion of:
    C.stagX_x = zeros (datos.Nx,datos.Ny); %X-staggered mesh for x-axis
    C.stagX_y = zeros (datos.Nx,datos.Ny); %X-staggered mesh for y-axis
    C.stagY_x = zeros (datos.Nx,datos.Ny); %Y-staggered mesh for x-axis
    C.stagY_y = zeros (datos.Nx,datos.Ny); %Y-staggered mesh for y-axis
    
    %Filling the meshes with each function
    [ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( datos );
    [ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( datos, C );
    [ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( datos, C );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %{
    figure; hold on;
    plot(C.Coll_x , C.Coll_y,'k.','MarkerSize', 35) ; 
    plot(C.stagX_x , C.stagX_y,'b.','MarkerSize',20) ;
    plot(C.stagY_x , C.stagY_y,'g.') ;
    ylim ([-0.1 1.1])
    xlim ([-0.1 1.1])
    %}    
    %Collocated Mesh
%{
    V1x = [C.dx(1) datos.L-C.dx(datos.Nx)];
    V1y = [C.dy(1) datos.H-C.dy(datos.Ny) ]; 
    V2x = [ datos.L-C.dx(datos.Nx)   datos.L-C.dx(datos.Nx) ]; 
    V2y = [ datos.H-C.dy(datos.Ny)   datos.H-C.dy(datos.Ny)];
    V3x = [C.dx(1) C.dx(1)];
    V3y = [C.dy(1) C.dy(1)];
   
    V5 = [C.dx(1)+C.dx(2) C.dx(1)+C.dx(2)];
    V6 = [C.dy(1) datos.H-C.dy(datos.Ny) ]; 
    V7 = [C.dx(1)+C.dx(2)+C.dx(3) C.dx(1)+C.dx(2)+C.dx(3)];
    
    V8 = [C.dx(1) datos.H-C.dy(datos.Ny) ];
    V9 = [C.dy(1)+C.dy(2) C.dy(1)+C.dy(2) ]; 
    V10 = [C.dy(1)+C.dy(2)+C.dy(3) C.dy(1)+C.dy(2)+C.dy(3)];
    %}
    %{
    if malla == 1
        hold on;
        plot(V1x,V2y,'k');  plot(V3x,V1y,'k'); plot(V2x,V1y,'k'); plot(V1x,V3y,'k');
        plot(V5,V6,'-.k'); plot(V7,V6,'-.k');
        plot(V8,V9,'-.k'); plot(V8,V10,'-.k');
        title('Collocated mesh'); 
    end
    
    % Staggered Mesh X
    if malla == 2
        V1x = V1x - C.dx(1)/2; V3x = V3x - C.dx(1)/2;
        V2x = V2x - C.dx(1)/2; 
    
        plot(V1x,V2y,'b');  plot(V3x,V1y,'b'); plot(V2x,V1y,'b'); plot(V1x,V3y,'b');
        plot(V5- C.dx(1)/2,V6,'-.b'); plot(V7- C.dx(1)/2,V6,'-.b');
        plot(V8- C.dx(1)/2,V9,'-.b'); plot(V8- C.dx(1)/2,V10,'-.b');
        title('Staggered mesh x-direction'); 
    end
        
    % Staggered Mesh Y
    if malla == 3
        V1y = V1y - C.dy(1)/2; V3y = V3y - C.dy(1)/2;
        V2y = V2y - C.dy(1)/2; 
        
        plot(V1x,V2y,'g');  plot(V3x,V1y,'g'); plot(V2x,V1y,'g'); plot(V1x,V3y,'g');
        plot(V5,V6- C.dy(1)/2,'-.g'); plot(V7,V6- C.dy(1)/2,'-.g');
        plot(V8,V9- C.dy(1)/2,'-.g'); plot(V8,V10- C.dy(1)/2,'-.g');
        title('Staggered mesh y-direction'); 
    end
    
    fprintf('---> Mallas calculada\n')      
   %}
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collocated_mesh function build the collocated mesh with no uniform meshes
% (computing delta X and delta Y for every type of distribution)
function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( datos )
    
    % Computing Delta X
    for i = 1:datos.Nx
            dx(i) = spacing_no_uniform( datos.gamma, datos.L, i, datos.Nx, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.L, i-1, datos.Nx, datos.uniform );
    end 
    
    % Computing Delta Y
    for j = 1:datos.Ny
            dy(j) = spacing_no_uniform( datos.gamma, datos.H, j, datos.Ny, datos.uniform ) - spacing_no_uniform( datos.gamma, datos.H, j-1, datos.Ny, datos.uniform );
    end 
    
    % Computing the X-position of COLLOCATED MESH
    X(1) = dx(1)/2; % X is inicialized in the average center of both VC faces
    for i = 2:datos.Nx
        X(i) = X(i-1) + dx(i-1); % Delta X is added in each defined position from the second one
    end
    
    
    % Computing the Y-position of COLLOCATED MESH
    Y(1) = dy(1)/2; % Y is inicialized in the average center of both VC faces
    for j = 2:datos.Ny
        Y(j) = Y(j-1) + dy(j-1); % Delta Y is added in each defined position from the second one
    end
    
    % Finally, the whole mesh is joined using as reference the coordenate
    % axis. 
    for i=1:datos.Nx
        for j=1:datos.Ny
            Coll_X(i,j) = X(j);
            Coll_Y(i,j) = Y(i);
        end
    end
    
 function div_u = divergencia_u(datos, u, v, nodal_mesh, num)
%This function computes the divergence of the velocity term in each CV.?
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y

Vx = datos.Vx;
Vy = datos.Vy;
Delta = datos.L/Vx;


 div_u = zeros(num,1); %This is vector ??u_p (divergence of velocity in each CV).

    for i=1:num
        
       [j k] = find(nodal_mesh == i);  %With this, one finds the coordinates of the CV.    
     
       u_p = u(k,j);
       v_p = v(k,j);
       
    if k-1 == 0 %Check if node is on the left edge
       u_w = u(j,Vx);
    else
       u_w = u(j,k-1); 
    end
    
    if j == Vx %check if node is on the lower edge 
     v_s = v(k,1);
    else
     v_s = v(k,j+1);
    end
    
    div_u(i) = Delta*(u_p-u_w+v_p-v_s);
    
      
end


% Stagg_mesh_x function build the staggered-x mesh with no uniform meshes
% (computing delta X for every type of distribution)
function [ X, Y ] = STAGG_MESH_X ( datos, C ) 
    
    % A x-direction vector is inicialized (with Nx size)
    vector_x = zeros( datos.Nx ,1 );
    
    % The mesh starts in the first position (0)
    vector_x(1) = 0;
    
    for i = 2:datos.Nx
        vector_x(i) = vector_x(i-1) + C.dx(i); 
        % From the first position, the incremental calculated before is 
        % added in each point. Adding the same delta that in the collocated
        % mesh assures the same VC length.
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %vector_x(datos.Nx+1) = datos.L;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % The X-Position is computed related to a vertical position. 
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = vector_x(i);
            Y(i,j) = C.Y(j);
       end
    end
   
    
    
    function [ X, Y ] = STAGG_MESH_Y ( datos, C )
    % This function works as the function STAGG_MESH_X
    
    vector_y = zeros( datos.Ny ,1 );
    
    vector_y(1) = 0; 
    
    for i = 2:datos.Ny
        vector_y(i) = vector_y(i-1) + C.dy(i);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %vector_y(datos.Ny) = datos.L;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    for i = 1:datos.Nx
       for j = 1:datos.Ny
            X(i,j) = C.X(i);
            Y(i,j) = vector_y(j);
       end
    end
    
   

% Spaging_no_uniform computes the Deltas (separations between VC faces)
% for every type of distribution    
function delta = spacing_no_uniform (gamma, L, i, N, uniform )

    % The reader should note that this function needs to know if the mesh
    % is uniform or not, creating a more efficient solution. 
    if uniform == true
        % L is the total mesh length. i is the position. N is the number of
        % nodes.
        delta = L * i / N;
    else
        delta = L*0.5*(1.0+tanh(gamma*(2.0*i/N-1.0))/tanh(gamma));
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

function du = diffusion_u(u,L)

N = size(u,2);
M = size(u,1);
dx = L/N;
dy = L/M;
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

for i = 2:N-1
    for j=2:M-1
        %Definim els valors de les velocitats als nodes W,E,S,N,P
        uE = u(i+1,j);
        uP = u(i,j);
        uW = u(i-1,j);
        uN = u(i,j+1);
        uS = u(i,j-1);

        du(i-1,j-1) = (uE - uP) / dx * Se + (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn + (uP-uS)/dy * Ss;
    end
end

function [p_gradX p_gradY] = gradient_p(datos, p, nodal_mesh, num)
%This function computes the divergence of the velocity term in each CV.?
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y

Vx = datos.Vx;
Vy = datos.Vy;
Delta_x = datos.L/Vx;
Delta_y = datos.L/Vy

%Transformem el vector de pseudo pressions en una matriu amb les pressions
%col·locades igual que en la matriu nodal mesh:
P_matrix = zeros(Vx,Vy);

for i=1:Vy
    for j=1:Vx
        P_matrix(i,j) = p(nodal_mesh(i,j));
    end
end

%Vectors to store the pressure gradient in x and y
p_gradX = zeros(Vx,Vy); 
p_gradY = zeros(Vx,Vy);

    for i=1:num
        
      [j k] = find(nodal_mesh == i);  %With this, one finds the coordinates of the CV. 
      
      p_p = P_matrix(j,k);
       
    if k == Vx %Check if node is on the right edge
       p_e = P_matrix(j,1);
    else
       p_e = P_matrix(j,k+1); 
    end
    
    if j == 1 %check if node is on the upper edge 
       p_n = P_matrix(k,Vy);
    else
       p_n = P_matrix(k,j-1);
    end
    
    p_gradX(j,k) = (p_e - p_p)/Delta_x;
    p_gradY(j,k) = (p_n - p_p)/Delta_y;
    
      
end
