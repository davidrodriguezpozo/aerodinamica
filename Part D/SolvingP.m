function [P P_prev dif_P] = SolvingP (ae, aw, an, as, ap, bp, P_prev, P, dif_P, d)
        
    delta = 1e-7;
    while dif_P>=delta 
         P_prev = P;
        
        for i = 2:d.Nrow-1
            for j = 2:d.Ncol-1
                P(i,j) = (ae(i,j)*P(i,j+1) + aw(i,j)*P(i,j-1) + an(i-1,j)*P(i-1,j) + as(i+1,j)*P(i+1,j) + bp(i,j))/ ap(i,j);      
            end
        end

         for i = 1:d.Nrow
            P(i,1) = P(i,2);
            P(i,d.Ncol) = P(i,d.Ncol-1);
         end
         
         for j = 1:d.Ncol
            P(d.Nrow,j) = P(d.Nrow-1,j);
            P(1,j) = P(2,j);
         end
         
         P(2,2) = 0;
         dif_P = 0;
         
         for i = 1:d.Nrow
             for j = 1:d.Ncol
                if abs(P(i,j)-P_prev(i,j))> d.dif_P
                    dif_P = abs(P(i,j)-P_prev(i,j));
                end
             end
         end
       
    end