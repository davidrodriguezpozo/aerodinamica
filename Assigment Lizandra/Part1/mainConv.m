% Main code for Assignament 2.
%   Barrachina, Victor
%   El Jarari, Younes
%   Royo, Enric

clear all; close all; 
clc;

run ('AirfoilData.m')

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

% Preprocess
N_w = b_w * 10;
%alpha_airfoil_w_vec = deg2rad([-3 0 3 6]);
alpha_airfoil_w_vec = deg2rad([3]);

Airfoils = [{'NACA_0009_N_16.txt'};{'NACA_0009_N_32.txt'};...
    {'NACA_0009_N_64.txt'};{'NACA_0009_N_128.txt'};...
    {'NACA_0009_N_256.txt'};{'NACA_0009_N_512.txt'}];

Ns = [16;32;64;128;256;512];

for p = 1:length(Airfoils)
    
X_airfoil_w = importdata (cell2mat(Airfoils(p)));
X_airfoil_w = X_airfoil_w (:,2:3);
N_airfoil_w = size(X_airfoil_w,1)-1;
c_w_vec = c_w1-linspace(0,N_w,N_w).*(c_w1-c_w2)/N_w;
theta_w_vec = 0-linspace(0,N_w,N_w).*(0-theta_w)/N_w;

%alpha_zl_w = compute_alphazl (N_airfoil_w, X_airfoil_w, rho, c_w)

    alpha_airfoil_w = alpha_airfoil_w_vec(1);
    alpha_ef_airfoil_w = alpha_airfoil_w;
    c_w = 1;

% Process
[cl_airfoil_w] = compute_cl (alpha_ef_airfoil_w, N_airfoil_w, X_airfoil_w, rho, c_w)
cl_vect (p) = cl_airfoil_w;

end

Airfoils = [{'NACA_2412_N_16.txt'};{'NACA_2412_N_32.txt'};...
    {'NACA_2412_N_64.txt'};{'NACA_2412_N_128.txt'};...
    {'NACA_2412_N_256.txt'};{'NACA_2412_N_512.txt'}];

for p = 1:length(Airfoils)
    
X_airfoil_w = importdata (cell2mat(Airfoils(p)));
X_airfoil_w = X_airfoil_w (:,2:3);
N_airfoil_w = size(X_airfoil_w,1)-1;
c_w_vec = c_w1-linspace(0,N_w,N_w).*(c_w1-c_w2)/N_w;
theta_w_vec = 0-linspace(0,N_w,N_w).*(0-theta_w)/N_w;

%alpha_zl_w = compute_alphazl (N_airfoil_w, X_airfoil_w, rho, c_w)

    alpha_airfoil_w = alpha_airfoil_w_vec(1);
    alpha_ef_airfoil_w = alpha_airfoil_w;
    c_w = 1;

% Process
[cl_airfoil_w] = compute_cl (alpha_ef_airfoil_w, N_airfoil_w, X_airfoil_w, rho, c_w)
cl_vecw (p) = cl_airfoil_w;

end



cl_final = cl_vect(end);
error_clt = abs(cl_vect-cl_final);
cl_final = cl_vecw(end);
error_clw = abs(cl_vecw-cl_final);

%figure;
%semilogx(1./Ns,cl_vec);

figure; 
loglog(1./Ns, error_clw, '-o'); hold on;
loglog(1./Ns, error_clt, '-o'); hold on;
%loglog(1./Ns, error_cdi, '-o');
%loglog(1./Ns, error_cdv, '-o');
grid minor; 
xlabel('$log_{10}$ $\frac{1}{N}$','Fontsize',16);
ylabel('$log_{10}$Error','Fontsize',16);
h=linspace(1/300, 0.1);
e=h.^2;
loglog(h,e,'r');
title('Logarithmic error of aerodynamic coefficients','Fontsize',18);

legend({'Main wing Lift coefficient $c_L$','Tail wing Lift coefficient $c_L$','$h^2$'},'Interpreter','latex','Fontsize',14,'Location','Northwest');




