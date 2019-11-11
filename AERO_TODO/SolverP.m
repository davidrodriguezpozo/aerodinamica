% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function [P P_prev dif_P] = SolverP (ae, aw, an, as, ap, bp, P_prev, P, dif_P, d)
        
    delta = 1e-7;
    dif_P = 1;
    while dif_P>=delta 
         P_prev = P;
        
        for i = 2:d.Nx-1
            for j = 2:d.Ny-1
                P(i,j) = (ae(i,j)*P(i,j+1) + aw(i,j)*P(i,j-1) +...
                    an(i,j)*P(i+1,j) + as(i,j)*P(i-1,j) + bp(i,j))/ ap(i,j);      
            end
        end

        P = haloupdate(P);
        dif_P = 0;
        
         
         for i = 1:d.Nx
             for j = 1:d.Ny
                if abs(P(i,j)-P_prev(i,j))> dif_P
                    dif_P = abs(P(i,j)-P_prev(i,j));
                end
             end
         end
       
    end