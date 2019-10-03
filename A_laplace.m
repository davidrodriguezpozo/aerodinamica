
Nx_vector = 3; 
Ny_vector = 3;

Vx = Nx_vector;

datos = INPUT(Vx,Vx);
C = meshes (datos, datos.malla);
A = A_matrix(Vx,Vx);


function A = A_matrix(Vx,Vy)
    nodal_mesh = zeros(Vx,Vy);
    num = Vx*Vy;
    ind=1;
    for i=Vx:-1:1
        for j=1:Vy
            nodal_mesh(i,j)=ind;
            ind = ind+1;
        end
    end
   
    A = zeros(num,num);
    
    for i=1:num
        %Diagonal term:
        A(i,i) = -4;
        
       [k j] = find(nodal_mesh == i);      
       %NORD:
       if k-1 == 0
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
           A(i, nodal_mesh(1,j)) = 1;
       else
           A(i, nodal_mesh(k,j+1)) = 1;
       end
    end
end


function datos = INPUT(Vx,Vy)

    datos.Nx = Vx+2;
    datos.Ny = Vy+2;
    
    datos.L = 1;
    datos.H = 1;
    
    datos.uniform = true;
    datos.gamma = 1;
    
    datos.malla = 2;
    
    datos.F = 1;
end