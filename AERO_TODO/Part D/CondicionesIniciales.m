function [S] = CondicionesIniciales (d, C)

fprintf('\n-( 2 )- Condiciones iniciales establecidas \n');

S.P = zeros(d.Nrow, d.Ncol);
S.P_prev = ones(d.Nrow, d.Ncol);

S.U = zeros(d.Nrow_sx, d.Ncol_sx);
S.U_pred = zeros(d.Nrow_sx, d.Ncol_sx);
S.U_prev = ones(d.Nrow_sx, d.Ncol_sx);
S.R_u = zeros(d.Nrow_sx, d.Ncol_sx);
S.R_u_prev = zeros(d.Nrow_sx, d.Ncol_sx);

S.V = zeros(d.Nrow_sy, d.Ncol_sy);
S.V_pred = zeros(d.Nrow_sy, d.Ncol_sy);
S.V_prev = ones(d.Nrow_sy, d.Ncol_sy);
S.R_v = zeros(d.Nrow_sy, d.Ncol_sy);
S.R_v_prev = zeros(d.Nrow_sy, d.Ncol_sy);

S.time = 0;
S.dt = d.dt_inicial;

for j = 1:d.Ncol_sx 
    U(1,j) = 1;
end
