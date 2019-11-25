function [r_1 r_2 theta_1 theta_2] = compute_r_theta (X,Xc,i,j);
        x_cp_i = Xc(i,1);   z_cp_i = Xc(i,2);
        x_j1 = X(j,1);      z_j1 = X(j,2);
        x_j2 = X(j+1);      z_j2 = X(j+1,2);

        r_1 = sqrt((x_j1-x_cp_i)^2 + (z_j1 - z_cp_i)^2);
        r_2 = sqrt((x_j2-x_cp_i)^2 + (z_j2 - z_cp_i)^2);
        theta_1 = asin((z_cp_i - z_j1) / (r_1));
        theta_2 = asin((z_cp_i - z_j2) / (r_2));
end