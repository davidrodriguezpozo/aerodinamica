function [c_l, c_l_alpha, c_m_14] = compute_coefficients (alpha_ef, N, X, rho, c)
    
    ca = cos(alpha_ef); sa = sin(alpha_ef);
    [Xc, n, t] = control_points (X,N);
    
%     plot (X(:,1),X(:,2))
%     hold on
% 
%     plot (Xc(:,1),Xc(:,2),'*')
%     plot (X(:,1),X(:,2),'o')

    U_inf = 1;
    %Q_inf = 1/2*U_inf^2*rho*[ca sa]; %cambio
    Q_inf = [ca sa];
    
    
    
    for i=1:N
        b(i,1) = -Q_inf * t(i,:)';
        for j=1:N
            % Compute the angle of the planes j
            x_j1 = X(j,1);      z_j1 = X(j,2);
            x_j2 = X(j+1);      z_j2 = X(j+1,2);
            [sa_j, ca_j, lj] = compute_angles(x_j1,x_j2,z_j1,z_j2);
            
            l(j,1) = lj;
            
            % Local coordinates of control points "i" in the local frame of
            % panel "j"
            x_ci_pan_j = (Xc(i,1)-X(j,1)) * ca_j - (Xc(i,2)-X(j,2)) * sa_j;
            z_ci_pan_j = (Xc(i,1)-X(j,1)) * sa_j + (Xc(i,2)-X(j,2)) * ca_j;
            
            % Induced velocity at "i"
            [r_1, r_2, theta_1, theta_2] = compute_r_theta (X,Xc,i,j);
            
            u_i_pan_j = (theta_2 - theta_1) / (2*pi);
            w_i_pan_j = log((r_2)^2 / (r_1)^2) / (4*pi);
            
            % Change the velocity back to global coordinates
            u_i = u_i_pan_j * ca_j + w_i_pan_j * sa_j;
            w_i = - u_i_pan_j * sa_j + w_i_pan_j * ca_j;
            V_i = [u_i; w_i];
            
            % Compute the influence ceofficients a_ij
            a(i,j) = [u_i w_i] * [ca_j -sa_j]';
        end
        C_p (i,1) = 1 - ((norm(V_i))/(norm(Q_inf)))^2;
        c_m0 (i,1) = 0;
    end
    for i=1:N
        a(i,i) = -.5; % CAMBIO a negativo
    end
    
    % Kutta condition
    i = 3;
    a(i,:) = 0; b(i,1) = 0;
    a(i,1) = 1; a(i,N) = 1;
    
    gamma = a\b;
    %gamma = inv(a)\b; 
    
    L = rho * Q_inf(1) * sum(gamma.*l);
    c_l =  L / (.5 * Q_inf(1)^2 * rho *c);
    %U_inf = 1; %Carlos
    %c_l = 2*sum(gamma.*l)/(U_inf*c); %Carlos
    
    c_l_alpha = 0;
    c_m_14 = 0;
    
end