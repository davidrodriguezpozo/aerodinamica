function S = Solve_NS (d, C, S, cp)

fprintf('\n-( 3 )- Calculando soluci�n... \n');

time = d.dt_inicial;
dt = d.dt_inicial;

time_steps = 0;

dx = d.L/d.Nx;

dif_V = d.dif_V;

while dif_V >= d.delta_V && time < d.tfinal 
    
    S.V_prev = S.V;
    %S.U_prev = S.U;
    
    % Update 
    
    S.R_u_prev = S.R_u; 
    S.R_v_prev = S.R_v; 
    
    [S.R_u S.R_v] = factor_R  (d, S.U, S.V);
    
    for i=2:d.Nrow_sx-1
        for j=2:d.Ncol_sx-1
            S.U_pred (i,j) = S.U(i,j) + dt*(1.5*S.R_u(i,j) - 0.5*S.R_u_prev(i,j))/(d.rho*dx*dx);
        end
    end
    
    for i=2:d.Nrow_sy-1
        for j=2:d.Ncol_sy-1
            S.V_pred (i,j) = S.V(i,j) + dt*(1.5*S.R_v(i,j) - 0.5*S.R_v_prev(i,j))/(d.rho*dx*dx);
        end
    end
    
    for i=2:d.Nrow-1
        for j=2:d.Ncol-1
            cp.bp (i,j) = -d.rho*dx*( S.U_pred(i-1,j) - S.U_pred(i-1,j-1) + S.V_pred(i-1,j-1) - S.V_pred(i,j-1))/dt;
        end
    end
    
    dif_P = 1;
    
    [S.P S.P_prev dif_P] = SolvingP (cp.ae, cp.aw, cp.an, cp.as, cp.ap, cp.bp, S.P_prev, S.P, dif_P, d);
    
    for i = 1:d.Nrow_sx
        for j = 1:d.Ncol_sx-1
            S.U(i,j) = S.U_pred(i,j) - dt*( S.P(i+1,j+1) - S.P(i+1,j))/(2*dx);
        end
    end
    
    for i = 2:d.Nrow_sy-1
        for j = 2:d.Ncol_sy-1
            S.V(i,j) = S.V_pred(i,j) - dt*( S.P(i,j+1) - S.P(i+1,j+1))/(2*dx);
        end
    end
    
    % BOUNDARY CONDITIONS
    
    S.U(d.Nrow_sx,:) = zeros(1,d.Ncol_sx);
    S.U(:,1) = zeros(d.Nrow_sx,1); 
    S.U(:,d.Ncol_sx) = zeros(d.Nrow_sx,1);
    
    S.V(1,:) = zeros(1,d.Ncol_sy); 
    S.V(d.Nrow_sy,:) = zeros(1,d.Ncol_sy);
    S.V(:,1) = zeros(d.Nrow_sy,1); 
    S.V(:,d.Ncol_sy) = zeros(d.Nrow_sy,1);
    
    S.U(1,:) = ones(1,d.Ncol_sx); 
    
    dt = CalculateTimeStep (d, S.U, S.V);
    time = time + dt;
    
    if time_steps >= 40
        fprintf('Simulation time: %f \n',time);
        time_steps = 0;
    end
    
    time_steps = time_steps + 1;
    dif_V = 0;
    
    for i = 1:d.Nrow_sx
        for j = 1:d.Ncol_sx
            if abs (S.U(i,j)-S.U_prev(i,j))> dif_V
                dif_V = abs(S.U(i,j)-S.U_prev(i,j));
            end
        end
    end
    
    
end

disp('----> Soluci�n calculada \n')
