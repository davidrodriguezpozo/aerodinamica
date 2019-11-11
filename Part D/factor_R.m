function [R_u R_v] = factor_R  (d, U, V)
    
    dx =d.L/d.Nx;
    
    R_u = zeros (d.Nrow_sx,d.Ncol_sx);
    R_v = zeros (d.Nrow_sy,d.Ncol_sy);
    
    for i=2:d.Nrow_sx-1
        for j = 2:d.Ncol_sx-1
            R_u(i,j)=d.mu*(U(i,j+1)+U(i,j-1)+U(i+1,j)+U(i-1,j)-4*U(i,j))- d.rho*0.25*dx*((U(i,j)+U(i,j+1))*(U(i,j)...
            +U(i,j+1))-(U(i,j)+U(i,j-1))*(U(i,j)+U(i,j-1))+(U(i,j)+U(i-1,j))*(V(i,j-1)+V(i,j))- (U(i,j)+U(i+1,j))*(V(i+1,j-1)+V(i+1,j)));
        end
    end
    
    %{
    for i=2:d.Nrow_sy-1
        for j = 2:d.Ncol_sy-1
            R_v(i,j)=d.mu*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j))-d.rho*0.25*dx*((U(i,j+1)+U(i-1,j+1))*(V(i,j)+V(i,j+1))-(U(i,j)+U(i-1,j))*(V(i,j)+V(i,j-1))+(V(i,j)+V(i-1,j))*(V(i,j)+V(i-1,j))- (V(i,j)+V(i+1,j))*(V(i,j)+V(i+1,j)));
        end
    end
    %}
    
    for i=2:d.Nrow_sy-1
        for j = 2:d.Ncol_sy-1
            R_v(i,j)=d.mu*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j))-d.rho*0.25*dx*((U(i,j+1)+U(i-1,j+1))*(V(i,j)+V(i,j+1))-(U(i,j)+U(i-1,j))*(V(i,j)+V(i,j-1))+(V(i,j)+V(i-1,j))*(V(i,j)+V(i-1,j))- (V(i,j)+V(i+1,j))*(V(i,j)+V(i+1,j)));
        end
    end
