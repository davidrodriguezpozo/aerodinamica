function d = INPUT

d.Re = 100;

d.Nx = 20;
d.Ny = 20;

d.Ncol = d.Nx + 2; % Numero de nodos en x-direction
d.Nrow = d.Ny + 2; % Numero de nodos en x-direction

d.Ncol_sx = d.Ncol - 1;
d.Nrow_sx = d.Nrow - 2;

d.Ncol_sy = d.Ncol - 2;
d.Nrow_sy = d.Nrow - 1;

d.L = 1;
d.H = 1;

d.tfinal = 50;

d.dt_inicial = 0.001;
d.rho = 1;
d.mu = 1/d.Re;

d.dif_V = 100;
d.dif_P = 100;

d.delta_P = 1e-7;
d.delta_V = 1e-7;
