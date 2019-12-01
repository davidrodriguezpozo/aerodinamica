function [c_l] = compute_cl (alpha_ef, N, X, rho, c)
    ca = cos(alpha_ef); sa = sin(alpha_ef);
    [Xc, n, t] = control_points (X,N);

    Q_inf = [ca sa];

    for i=1:N
        b(i,1) = -Q_inf * t(i,:)';
        for j=1:N
            % Compute the angle of the planes j
            x_j1 = X(j,1);      z_j1 = X(j,2);
            x_j2 = X(j+1,1);    z_j2 = X(j+1,2);
            [sa_j, ca_j, lj] = compute_angles(x_j1,x_j2,z_j1,z_j2);
            
            l(j,1) = lj;
            
            % Local coordinates of control points "i" in the local frame of
            % panel "j"
            x_ci_pan_j = (Xc(i,1)-X(j,1))*ca_j - (Xc(i,2)-X(j,2))*sa_j;
            z_ci_pan_j = (Xc(i,1)-X(j,1))*sa_j + (Xc(i,2)-X(j,2))*ca_j;
            
            % Induced velocity at "i"
           
            r22_r12 = ((x_ci_pan_j-lj)^2+z_ci_pan_j^2)/(x_ci_pan_j^2+z_ci_pan_j^2);
            theta_1 = atan2(z_ci_pan_j,x_ci_pan_j);
            theta_2 = atan2(z_ci_pan_j,(x_ci_pan_j-lj));
        
            u_i_pan_j = (theta_2 - theta_1) / (2*pi);
            w_i_pan_j = log(r22_r12) / (4*pi);
            
            % Change the velocity back to global coordinates
            u_i = u_i_pan_j * ca_j + w_i_pan_j * sa_j;
            w_i =-u_i_pan_j * sa_j + w_i_pan_j * ca_j;
            V_i = [u_i; w_i];
            
            % Compute the influence ceofficients a_ij
            a(i,j) = dot([u_i, w_i], t(i,:));
        end
        a(i,i) = -.5;
    end
     
    % Kutta condition
    i = N/4;
    a(i,:) = 0; b(i,1) = 0;
    a(i,1) = 1; a(i,N) = 1;
    
    gamma = a\b;
    
    L = rho * Q_inf(1) * sum(gamma.*l);
    c_l = L / (.5 * Q_inf(1)^2 * rho *c);
    
end