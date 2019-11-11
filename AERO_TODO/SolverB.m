function [u, v, P, time, errorP] = SolverB(datos, C)

    [u v R_u R_v u_p v_p u_prev v_prev] = CondicionesIniciales(datos, C);
    
    P = zeros(datos.Nx);
    P_prev = zeros(datos.Nx);
    delta_T = 0.001;
    
    [cp] = Coeff_Pressure (datos);
    
    Nx = datos.Nx;
    Ny = datos.Ny;
    
    time = 0; final_time = 2;
    dif_U = 1;
    dif_V = 1;
    dif_P = 1;
    delta_V = 1e-7;
    
    while time <= final_time && dif_U > delta_V && dif_V > delta_V
    
    R_u_prev = R_u;
    R_v_prev = R_v;
    
    u_prev = u;
    v_prev = v;
    
    %datos.F = exp(-4*pi*datos.mu/datos.rho*time);
    
    [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v);
    
    %Calculem R: R = - Conv(u)/V + Diff(u)/V
    Vol = (datos.L/datos.Vx)*(datos.H/datos.Vy);
    R_u = -(conv_u/Vol) + datos.mu*(diff_u/Vol);
    R_v = -(conv_v/Vol) + datos.mu*(diff_v/Vol);
    
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1
        u_p(i,j) = u(i,j) + delta_T*(3/2*R_u(i,j) - 1/2*R_u_prev(i,j));
        v_p(i,j) = v(i,j) + delta_T*(3/2*R_v(i,j) - 1/2*R_v_prev(i,j));
        end
    end
    
    u_p = haloupdate(u_p);
    v_p = haloupdate(v_p);
    
    dx = datos.L/datos.Vx;
    dy = datos.H/datos.Vy;
    rho = datos.rho;
    mu = datos.mu;

    for i=2:Nx-1
        for j=2:Ny-1
            cp.bp(i,j)=-1/delta_T*rho*(u_p(i,j)*(dy)-u_p(i,j-1)*(dy)+v_p(i,j)*(dx)-v_p(i-1,j)*(dx));
        end
    end

    [P P_prev dif_P] = SolverP (cp.ae, cp.aw, cp.an, cp.as, cp.ap, cp.bp, P_prev, P, dif_P, datos);
    
    P = haloupdate(P);
    
    p_analytic = Analytic_Pressure (datos, C, time);

    for i=2:Nx-1
        for j=2:Nx-1
 
        u(i,j)=u_p(i,j)-delta_T/rho*(P(i,j+1)-P(i,j))/(dx);
        v(i,j)=v_p(i,j)-delta_T/rho*(P(i+1,j)-P(i,j))/(dy);
 
        end
    end
    
    u = haloupdate(u);
    v = haloupdate(v);
    
    dif_U = 0;
    dif_V = 0;
    
    for i = 1:Nx
        for j = 1:Ny
            if abs(u(i,j)-u_prev(i,j))> dif_U
                dif_U = abs(u(i,j)-u_prev(i,j));
            end
            
            if abs(v(i,j)-v_prev(i,j))> dif_V
                dif_V = abs(v(i,j)-v_prev(i,j));
            end
        end
    end
    
    delta_T = TimeStep (datos,u,v);
    
    time = time + delta_T;
    
    datos.F = exp(-8*pi^2*datos.mu*time);

end

    errorP = ERROR(p_analytic, P);