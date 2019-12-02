function plots(Nx_vector, error_du, error_dv, error_cu, error_cv)
quad = (1./Nx_vector).^2;


figure;
loglog(1./Nx_vector,quad); hold on;
%Provo ap lotejar tots els errors junts a ver que pasa
loglog(1./Nx_vector,0.5*(error_dv+error_du),'ro-'); hold on;
loglog(1./Nx_vector,0.5*(error_cv+error_cu),'ko-');
% xlim([xlim(1)-1 xlim(2)+1]);
% ylim([ylim(1)-1 ylim(2)+1]);
% xlim([1/Nx_vector(1)-0.1 1/Nx_vector(length(Nx_vector))+1]);
% ylim([min(quad)-0.5 max(error_du)+1]);
legend('$h^2$','Diffusive term','Convective term','Interpreter','latex');
title('Diffusive and convective terms error vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error');

figure;
loglog(1./Nx_vector,error_du,'o-'); hold on;
title('Diffusive term error (x-direction) vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error diffusive term (x-dir)');


% figure;
% loglog(1./Nx_vector,error_dv,'o-');
% title('Diffusive term error (y-direction) vs. mesh size','Fontsize',16);
% xlabel('Size of the mesh h');
% ylabel('Error diffusive term (y-dir)');
% 
% figure;
% loglog(1./Nx_vector,error_cu,'o-');
% title('Convective term error (x-direction) vs. mesh size','Fontsize',16);
% xlabel('Size of the mesh h');
% ylabel('Error convective term (x-dir)');
% 
% figure;
% loglog(1./Nx_vector,error_cv,'o-');
% title('Convective term error (y-direction) vs. mesh size','Fontsize',16);
% xlabel('Size of the mesh h');
% ylabel('Error convective term (y-dir)');
% 
% figure;
% loglog(1./Nx_vector,error_cu,'o-b'); hold on;
% loglog(1./Nx_vector,error_cv,'o-r');
% title('Convective term error vs. mesh size','Fontsize',18);
% xlabel('Size of the mesh h','Fontsize',16);
% ylabel('Error convective term','Fontsize',16);
% legend('x-direction','y-direction','Fontsize',16,'Location','southeast');
% 
% figure;
% loglog(1./Nx_vector,error_du,'o-b','Linewidth',1.5); hold on;
% loglog(1./Nx_vector,error_dv,'o-r');
% title('Diffusive term error vs. mesh size','Fontsize',18);
% xlabel('Size of the mesh h','Fontsize',16);
% ylabel('Error diffusive term','Fontsize',16);
% legend('x-direction','y-direction','Fontsize',16,'Location','southeast');

