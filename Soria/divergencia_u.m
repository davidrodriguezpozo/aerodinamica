function div_u = divergencia_u(datos, u, v, nodal_mesh, num)
%This function computes the divergence of the velocity term in each CV.�
% u : velocity values of staggered mesh-X amb HALO 
% v : velocity values of staggered mesh-Y

Vx = datos.Vx;
Vy = datos.Vy;
Delta = datos.L/Vx;


 div_u = zeros(num,1); %This is vector ?�u_p (divergence of velocity in each CV).

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


end


