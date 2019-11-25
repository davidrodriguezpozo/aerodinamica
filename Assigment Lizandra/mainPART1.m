clear all
close all
clc

run ('AirfoilData.m')

alpha_w_vec = deg2rad([-3 0 3 6]);
divs = [16 32 64 128 256 512];

%X_w = importdata ('NACA_2412_N_512.txt');

alpha_w = deg2rad(3);

for i=1:length(divs)
    for j=1:length(alpha_w_vec)

        X_w = Geometry(divs(i),'wing');
        alpha_w = alpha_w_vec(j);

        X_w = X_w (:,2:3);
        N_w = size(X_w,1)-1;

        alpha_ef_w = alpha_w;
        c_w = 1;

        % Process
        [cl_w, cl_alpha_w, c_m_14_w] = ...
        compute_coefficients (alpha_ef_w, N_w, X_w, rho, c_w);
        cl_w_vec (i,j) = cl_w;
        c_m_14_w_vec (i,j) = c_m_14_w;
        
    end
    
    for k = 1:length(alpha_w_vec)-1
        cl_alpha_w_vec(k) = (cl_w_vec (i,k+1) - cl_w_vec (i,k))/(alpha_w_vec(k+1) - alpha_w_vec(k));
    end

    cl_alpha_w(i) = sum(cl_alpha_w_vec)/length(cl_alpha_w_vec);
    
end

%figure;
%semilogx (divs,cl_w_vec);

figure;
error1 = abs(cl_w_vec(:,1)-cl_w_vec(end,1));
error2 = abs(cl_w_vec(:,2)-cl_w_vec(end,2));
error3 = abs(cl_w_vec(:,3)-cl_w_vec(end,3));
error4 = abs(cl_w_vec(:,4)-cl_w_vec(end,4));
semilogx (divs,error1); hold on;
semilogx (divs,error2); 
semilogx (divs,error3); 
semilogx (divs,error4); 
%plot (divs,error1); hold on;
%plot (divs,error2); 
%plot (divs,error3); 
%plot (divs,error4); 
ylim([min(error1)-0.05 max(error1)+0.05])

title('Error c_l for different AoA');

%% Per altres plots

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
   cl_alpha_w_vec(i) = (cl_w_vec (i+1) - cl_w_vec (i))/(alpha_w_vec(i+1) - alpha_w_vec(i));
end

cl_alpha_w = sum(cl_alpha_w_vec)/length(cl_alpha_w_vec);

figure;
plot (rad2deg(alpha_w_vec),alpha_w_vec,'-*')




