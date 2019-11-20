function [p_gradX p_gradY] = gradient_p(datos, p, nodal_mesh, num)
%This function computes the divergence of the velocity term in each CV.?
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y

Vx = datos.Nx;
Vy = datos.Ny;
Delta_x = datos.L/datos.Vx;
Delta_y = datos.L/datos.Vy;

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

    for i=2:num-1
        
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