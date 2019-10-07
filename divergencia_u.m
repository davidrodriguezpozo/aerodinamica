function u_p = divergencia_u(datos, C, u, v)
%This function computes the divergence of the velocity term in each CV.ç
% u : velocity values of staggered mesh-X
% v : velocity values of staggered mesh-Y

Vx = datos.Vx;
Vy = datos.Vy;

[nodal_mesh num] = nodalmesh(Vx,Vy);


 u_p = zeros(num,1); %This is vector ?·u_p (divergence of velocity in each CV).

    for i=1:num
        
       [k j] = find(nodal_mesh == i);  %With this, one finds the coordinates of the CV.    
     
       u_p = u(k,j);
       v_p = v(k,j);
       
    if k-1 == 0 %Check if node is on the left edge
       u_w = u(Vx,j);
    else
       u_w = u(k-1,j); 
    end
    
    if j-1 == 0 %check if node is on the upper edge 
     v_s = v(k,Vy);
    else
     v_s = v(k,j-1);
    end
    
    u_p(i) = C.dx(1)*(u_p-u_w+v_p-v_s);
    
      
end


end


