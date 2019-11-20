function plots(Nx_vector, error_du, error_dv, error_cu, error_cv)

figure;
loglog(1./Nx_vector,error_du,'o-');
title('Diffusive term error (x-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error diffusive term (x-dir)');

figure;
loglog(1./Nx_vector,error_dv,'o-');
title('Diffusive term error (y-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error diffusive term (y-dir)');

figure;
loglog(1./Nx_vector,error_cu,'o-');
title('Convective term error (x-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error convective term (x-dir)');

figure;
loglog(1./Nx_vector,error_cv,'o-');
title('Convective term error (y-direction) vs. mesh size');
xlabel('Size of the mesh h');
ylabel('Error convective term (y-dir)');