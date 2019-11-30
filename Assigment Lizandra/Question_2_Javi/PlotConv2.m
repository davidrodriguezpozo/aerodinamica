clc; clear all; close all;

format long
%% INPUT DATA
run('InputDataA2');

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');


disp('NO GROUND EFFECT RESULTS:')
Ground_effect=0;

Ns = [8 16 32 64 128 256 512];
M = length(Ns);

for i = 1:M
    N = Ns(i);
    [cl_vec(i) cdv_vec(i) cdi_vec(i) cd_vec(i)] = Task21(Ground_effect, h, N, alpha, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3);
end

cl_final = cl_vec(end);
error_cl = abs(cl_vec-cl_final);
cdv_final = cdv_vec(end);
error_cdv = abs(cdv_vec-cdv_final);
cdi_final = cdi_vec(end);
error_cdi = abs(cdi_vec-cdi_final);
cd_final = cd_vec(end);
error_cd = abs(cd_vec-cd_final);

%figure;
%semilogx(1./Ns,cl_vec);

figure; 
loglog(1./Ns, error_cl, '-o'); hold on;
loglog(1./Ns, error_cdi, '-o');
loglog(1./Ns, error_cdv, '-o');
grid minor; 
xlabel('$log_{10}$ $\frac{1}{N}$','Fontsize',16);
ylabel('$log_{10}$Error','Fontsize',16);
h=linspace(1/1000, 0.15);
e=h.^2;
loglog(h,e,'r');
title('Logarithmic error of aerodynamic coefficients','Fontsize',18);

legend({'Lift coefficient $C_L$','Induced Drag $C_{D_i}$','Viscous Drag $C_{D_v}$','$h^2$'},'Interpreter','latex','Fontsize',14,'Location','Northwest');




