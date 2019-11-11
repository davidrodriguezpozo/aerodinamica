function plotsA(Nx_vector, error_du, error_dv, error_cu, error_cv)

figure; 
loglog(1./Nx_vector,error_du,'o-');
title('Diffusive term error (x-direction) vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error diffusive term (x-dir)');

figure;
loglog(1./Nx_vector,error_dv,'o-');
title('Diffusive term error (y-direction) vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error diffusive term (y-dir)');

figure;
loglog(1./Nx_vector,error_cu,'o-');
title('Convective term error (x-direction) vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error convective term (x-dir)');

figure;
loglog(1./Nx_vector,error_cv,'o-');
title('Convective term error (y-direction) vs. mesh size','Fontsize',16);
xlabel('Size of the mesh h');
ylabel('Error convective term (y-dir)');

figure;
loglog(1./Nx_vector,error_cu,'o-b'); hold on;
loglog(1./Nx_vector,error_cv,'o-r');
title('Convective term error vs. mesh size','Fontsize',18);
xlabel('Size of the mesh h','Fontsize',16);
ylabel('Error convective term','Fontsize',16);
legend('x-direction','y-direction','Fontsize',16,'Location','southeast');

figure;
loglog(1./Nx_vector,error_du,'o-b','Linewidth',1.5); hold on;
loglog(1./Nx_vector,error_dv,'o-r');
title('Diffusive term error vs. mesh size','Fontsize',18);
xlabel('Size of the mesh h','Fontsize',16);
ylabel('Error diffusive term','Fontsize',16);
legend('x-direction','y-direction','Fontsize',16,'Location','southeast');