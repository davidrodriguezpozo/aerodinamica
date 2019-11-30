clc; clear all; close all;

format long
%% INPUT DATA
run('InputDataA2');

%set(groot, 'DefaultTextInterpreter','latex');
%set(groot, 'Defaultaxesticklabelinterpreter','latex');


disp('NO GROUND EFFECT RESULTS:')
Ground_effect=0;

Ns = [16 32 64 128 256 512];
M = length(Ns);

for i = 1:M
    N = Ns(i);
    [cl_vec(i) cdv_vec(i) cdi_vec(i) cd_vec(i)] = Task21(Ground_effect, h, N, alpha, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3);
end

cl_final = cl_vec(end);
error_cl = abs(cl_vec-cl_final);
cd_final = cd_vec(end);
error_cd = abs(cd_vec-cd_final);

a = 1:20; a= 0.1*a;
b = 1:20;
figure;
loglog(a,b);

%figure;
%semilogx(1./Ns,cl_vec);

figure; hold on;
loglog(1./Ns, error_cl, '-o');
%loglog(1./Ns, error_cd, '-o');
grid minor; 
xlabel('$log_{10}$ $\frac{1}{N}$','Interpreter','latex','Fontsize',14);
ylabel('$log_{10}$Error','Interpreter','latex','Fontsize',14);
h=linspace(1/1000, 0.1);
e=h.^2;
loglog(h,e,'r');
title('Logarithmic error of $C_{L}$','Interpreter','latex','Fontsize',20);

legend('Lift coefficient','h^2');




