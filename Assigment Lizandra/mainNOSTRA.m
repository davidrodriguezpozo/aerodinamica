function main
clear all
close all
clc

run ('AirfoilData.m')

alpha_w_vec = deg2rad([-3 0 3 6]);
divs = [16 32 64 128 256 512];


X = Geometry(divs(6),'wing');
alpha_w = alpha_w_vec(3);

X = X (:,2:3);
N = size(X,1)-1;

alpha_ef = alpha_w;
c = 1;

% Process
[cl, cl_alpha, c_m_14] = Coefficients (alpha_ef, N, X, rho, c);

end


function [cl, cl_alpha, cm_14] = Coefficients (alpha_ef, N, X, rho, c)
    ca = cos(alpha_ef); sa = sin(alpha_ef);
    [Xc, n, t] = control_points (X,N);
    
    for i=1:N
        x1 = X(i,1); x2 = X(i+1,1);
        z1 = X(i,2); z2 = X(i+1,2);
        li = sqrt((x1-x2)^2+(z1-z2)^2);
        
        %[sa ca] = compute_angles (x1,x2,z1,z2);
        %sa = (z2-z1) / (li); 
        % Cambio sa
        si = (z1-z2) / (li); 
        ci = (x2-x1) / (li);
        
        %Xc(i,1) = x1 + .5*li*ca; % Mid point
        %Xc(i,2) = z1 + .5*li*sa; % Mid point
        
        Xc(i,1) = (x2 + x1)/2; % Cambio
        Xc(i,2) = (z2 + z1)/2; % Cambio
        
        n(i,1) = si; n(i,2) = ci;
        t(i,1) = ci; t(i,2) = -si;
    end
    
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
            %[sa_j, ca_j, lj] = compute_angles(x_j1,x_j2,z_j1,z_j2);
            
            x1 = X(i,1); x2 = X(i+1,1);
            z1 = X(i,2); z2 = X(i+1,2);
            
            sa_j = (z1-z2) / (li); 
            ca_j = (x2-x1) / (li);

            lj = sqrt((x1-x2)^2+(z1-z2)^2);
            
            l(j,1) = lj;
            
            % Local coordinates of control points "i" in the local frame of
            % panel "j"
            x_ci_pan_j = (Xc(i,1)-X(j,1)) * ca_j - (Xc(i,2)-X(j,2)) * sa_j;
            z_ci_pan_j = (Xc(i,1)-X(j,1)) * sa_j + (Xc(i,2)-X(j,2)) * ca_j;
            
            % Induced velocity at "i"
            %[r_1, r_2, theta_1, theta_2] = compute_r_theta (X,Xc,i,j);
            
            x_cp_i = Xc(i,1);   z_cp_i = Xc(i,2);
            x_j1 = X(j,1);      z_j1 = X(j,2);
            x_j2 = X(j+1,1);      z_j2 = X(j+1,2);

            r_1 = sqrt((x_j1-x_cp_i)^2 + (z_j1 - z_cp_i)^2);
            r_2 = sqrt((x_j2-x_cp_i)^2 + (z_j2 - z_cp_i)^2);
            %theta_1 = asin((z_cp_i - z_j1) / (r_1));
            %theta_2 = asin((z_cp_i - z_j2) / (r_2));

            theta_1 = atan((z_j1 - z_cp_i) / (x_j1 - x_cp_i)); %Carlos
            theta_2 = atan((z_j2 - z_cp_i) / (x_j2 - x_cp_i)); %Carlos
            
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
    cl = L / (.5 * Q_inf(1)^2 * rho *c);
    %U_inf = 1; %Carlos
    c_l = 2*sum(gamma.*l)/(U_inf*c); %Carlos
    
    cl_alpha = 0;
    cm_14 = 0;
    
end



