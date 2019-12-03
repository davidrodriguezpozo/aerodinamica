function div_u = divergencia_u(datos, u, v)
%This function computes the divergence of the velocity term in each CV.�
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y amb HALO

% Li entrem la matriu amb el halo
Vx = datos.Vx;
Vy = datos.Vy;
Delta = datos.L/Vx;


%  div_u = zeros(num,1);
 div_u = zeros(Vx*Vy,1); %This is vector Delta_u_p (divergence of velocity in each CV).

%     for i=1:num
%         
%        [j k] = find(nodal_mesh == i);  %With this, one finds the coordinates of the CV.    
%      
%        u_p = u(k,j);
%        v_p = v(k,j);
%        
%     if k-1 == 0 %Check if node is on the left edge
%        u_w = u(j,Vx);
%     else
%        u_w = u(j,k-1); 
%     end
%     
%     if j == Vx %check if node is on the lower edge 
%      v_s = v(k,1);
%     else
%      v_s = v(k,j+1);
%     end
    for i=2:Vx+1
        for j=2:Vy+1
           u_p = u(i,j);
           v_p = v(i,j);
           u_w = u(i-1,j);
           v_s = u(i,j-1);
           
        div_u(i) = Delta*(u_p-u_w+v_p-v_s);
        end
      
    end


end


