function [r22_r12, theta_1, theta_2] = compute_r_theta (x,z,l);
        r22_r12 = ((x-l)^2+z^2)/(x^2+z^2);
        theta_1 = atan(z/x);
        theta_2 = atan(z/(x-l));
        
end