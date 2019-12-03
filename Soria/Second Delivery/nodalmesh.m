function [nodal_mesh, num] = nodalmesh(Vx,Vy)
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
end
