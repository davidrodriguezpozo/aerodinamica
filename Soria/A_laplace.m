function A = A_laplace(datos,C)

Vx = datos.Vx;
Vy = datos.Vy;

[nodal_mesh, num] = nodalmesh(Vx,Vy);


   
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
end

