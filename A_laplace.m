function A = A_laplace(datos,C)


Vx = datos.Vx;

A = A_matrix(Vx,Vx);

divA = divergence (A);


function divA = divergence (A)

divA = 0;

end



function A = A_matrix(Vx,Vy)
    nodal_mesh = zeros(Vx,Vy); %Matriu de tots els nodes dels VC.
    num = Vx*Vy; %Per sabre fins a quin index arribem
    ind=1;
    for i=Vx:-1:1
        for j=1:Vy
            nodal_mesh(i,j)=ind; %S'ompla la matriu de nodes amb els index necessaris;
            ind = ind+1;
        end
    end
   
    A = zeros(num,num); %Ser� la matriu que es multiplicar� per p.
    
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
    
    A(9,9) =  -4;
end


end

