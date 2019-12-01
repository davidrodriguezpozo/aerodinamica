% Main code for Assignament 2.
%   Barrachina, Victor
%   El Jarari, Younes
%   Royo, Enric

clear all; close all; 
clc;

run ('AirfoilData.m')

% Preprocess
N_w = b_w * 10;
alpha_airfoil_w_vec = deg2rad([-3 0 3 6]);
% alpha_airfoil_w_vec = deg2rad([3]);

X_airfoil_w = importdata ('NACA_0009_N_512.txt');
X_airfoil_w = X_airfoil_w (:,2:3);
N_airfoil_w = size(X_airfoil_w,1)-1;
c_w_vec = c_w1-linspace(0,N_w,N_w).*(c_w1-c_w2)/N_w;
theta_w_vec = 0-linspace(0,N_w,N_w).*(0-theta_w)/N_w;

%alpha_zl_w = compute_alphazl (N_airfoil_w, X_airfoil_w, rho, c_w)

for i = 1:length(alpha_airfoil_w_vec)
    alpha_airfoil_w = alpha_airfoil_w_vec(i);
    alpha_ef_airfoil_w = alpha_airfoil_w;
    c_w = 1;

% Process
[cl_airfoil_w, Vi, Cp, cmo] = compute_cl (alpha_ef_airfoil_w, N_airfoil_w, X_airfoil_w, rho, c_w)
cl_airfoil_w_vec (i) = cl_airfoil_w;

end

for i = 1:length(alpha_airfoil_w_vec)-1
   cl_alpha_w_vec(i) = (cl_airfoil_w_vec (i+1) - cl_airfoil_w_vec (i))/(alpha_airfoil_w_vec(i+1) - alpha_airfoil_w_vec(i));
end
cl_alpha_w = mean (cl_alpha_w_vec)

plot (rad2deg(alpha_airfoil_w_vec),cl_airfoil_w_vec,'-*')


    


