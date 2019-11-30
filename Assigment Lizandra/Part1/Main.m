% Main code for Assignament 2.
%   Barrachina, Victor
%   El Jarari, Younes
%   Royo, Enric

%clear all; close all; 
clc;

run ('AirfoilData.m')

% Preprocess
alpha_w_vec = deg2rad([-3 0 3 6]);

X_w = importdata ('NACA_2412_N_512.txt');
X_w = X_w (:,2:3);
N_w = size(X_w,1)-1;

for i = 1:length(alpha_w_vec)
    alpha_w = alpha_w_vec(i);
    alpha_ef_w = alpha_w;
    c_w = 1;

% Process
[cl_w, cl_alpha_w, c_m_14_w] = ...
compute_coefficients (alpha_ef_w, N_w, X_w, rho, c_w)
cl_w_vec (i) = cl_w;
end

for i = 1:length(alpha_w_vec)-1
   cl_alpha_w_vec = (cl_w_vec (i+1) - cl_w_vec (i))/(alpha_w_vec(i+1) - alpha_w_vec(i))
end

plot (rad2deg(alpha_w_vec),alpha_w_vec,'-*')

    

