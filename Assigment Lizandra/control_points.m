function [Xc n t] = control_points (X,N)
    for i=1:N
        x1 = X(i,1); x2 = X(i+1,1);
        z1 = X(i,2); z2 = X(i+1,2);
        li = sqrt((x1-x2)^2+(z1-z2)^2);
        
        %[sa ca] = compute_angles (x1,x2,z1,z2);
        %sa = (z2-z1) / (li); 
        % Cambio sa
        sa = (z1-z2) / (li); 
        ca = (x2-x1) / (li);
        
        Xc(i,1) = x1 + .5*li*ca;
        Xc(i,2) = z1 + .5*li*sa;
        
        n(i,1) = sa; n(i,2) = ca;
        t(i,1) = ca; t(i,2) = -sa;
    end
end