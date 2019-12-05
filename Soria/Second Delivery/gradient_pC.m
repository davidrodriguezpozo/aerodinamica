function [p_gradX, p_gradY, P_matrix_halo] = gradient_pC(datos, p, nodal_mesh,delta_t)
%This function computes the gradient of the pressure in each CV.
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y

Vx = datos.Vx;
Vy = datos.Vy;
Delta_x = datos.L/Vx;
Delta_y = datos.L/Vy;

%Transformem el vector de pseudo pressions en una matriu amb les pressions
%col·locades igual que en la matriu nodal mesh:
P_matrix = zeros(Vx,Vy);

for i=1:Vy
    for j=1:Vx
        k = nodal_mesh(i,j);
        P_matrix(i,j) = p(k); % P_matrix =  [Vx,Vy]
    end
end

P_matrix_halo = zeros(Vx+2,Vy+2);

for i=2:Vx+1 %Redimensionem el vector de pressions a [Nx,Ny]
    for j=2:Vy+1
    
    P_matrix_halo(i,j) = P_matrix(i-1,j-1);

    end
end

P_matrix_halo = haloupdate(P_matrix_halo); %Fem el halo de la P_matrix


%Vectors to store the pressure gradient in x and y
p_gradX = zeros(Vx+2,Vy+2); % [Nx x Ny]
p_gradY = zeros(Vx+2,Vy+2); % [Nx x Ny]

for i=2:Vx+1
    for j=2:Vy+1

% 
% %     for i=1:num
% %         
%       [j k] = find(nodal_mesh == i);  %With this, one finds the coordinates of the CV. 
%       
%       p_p = P_matrix(j,k);
%        
%     if k == Vx %Check if node is on the right edge
%        p_e = P_matrix(j,1);
%     else
%        p_e = P_matrix(j,k+1); 
%     end
%     
%     if j == 1 %check if node is on the upper edge 
%        p_n = P_matrix(k,Vy);
%     else
%        p_n = P_matrix(k,j-1);
%     end
    p_e = P_matrix_halo(i,j+1);
    p_p = P_matrix_halo(i,j);
    p_n = P_matrix_halo(i,j-1);
    
    p_gradX(i,j) = delta_t/datos.rho*(p_e - p_p)/Delta_x;
    p_gradY(i,j) = delta_t/datos.rho*(p_n - p_p)/Delta_y;

    end
end


end