function main_B

    set(groot, 'DefaultTextInterpreter','latex');
    set(groot, 'Defaultaxesticklabelinterpreter','latex');

    clear all 
    clc
    close all 


    Vx = 30; 
    Vy = 30;
    Re = 1000;

    datos = INPUT(Vx,Vy,Re);
    C = meshes (datos, datos.malla);

    datos.F = 1;

    u = zeros(datos.Nx, datos.Ny);
    v = zeros(datos.Nx, datos.Ny);
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1

            x1 = C.stagX_x(i,j);
            y1 = C.stagX_y(i,j);

            u(i,j) = datos.F*cos(2*pi*x1)*sin(2*pi*y1);

            x2 = C.stagY_x(i,j);
            y2 = C.stagY_y(i,j);

            v(i,j) = -datos.F*cos(2*pi*y2)*sin(2*pi*x2);

        end
    end

u = haloupdate(u);
v = haloupdate(v);

R_u = zeros(datos.Nx);
R_v = zeros(datos.Nx);
p = zeros(Vx^2,1);
p_prev = zeros(Vx^2,1);

conv_u = zeros(datos.Nx); 
diff_u = zeros(datos.Nx);
conv_v = zeros(datos.Nx);
diff_v = zeros(datos.Nx);

u_p = zeros(datos.Nx);
v_p = zeros(datos.Nx);
u_prev = zeros(datos.Nx);
v_prev = zeros(datos.Nx);

    ae = zeros(datos.Nx);
    aw = zeros(datos.Nx);
    an = zeros(datos.Nx);
    as = zeros(datos.Nx);
    ap = zeros(datos.Nx);
    bp = zeros(datos.Nx);

    P = zeros(datos.Nx);
    P_prev = zeros(datos.Nx);
    dif_P = 1;
    delta_T = 0.01;

    Nx = datos.Nx;
    Ny = datos.Ny;

    R_u_prev = R_u;
    R_v_prev = R_v;
    
    u_prev = u;
    v_prev = v;
    
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
    
    rho = 1;
    
    for i=2:Nx-1
        for j=2:Ny-1
            ae(i,j) = dy/dx;
            aw(i,j) = dy/dx;
            an(i,j) = dx/dy;
            as(i,j) = dx/dy; 
            bp(i,j)=-1/delta_T*rho*(u_p(i,j)*(dy)-u_p(i,j-1)*(dy)+v_p(i,j)*(dx)-v_p(i-1,j)*(dx));
            
        end
    end
    
    ap = ae + aw + an + as;
   
    [P P_prev dif_P] = SolverP (ae, aw, an, as, ap, bp, P_prev, P, dif_P, datos);
    
    P = haloupdate(P);
    time = delta_T;
    p_analytic = Analytic_Pressure (datos, C, time);
    
    rho = datos.rho;
    mu = datos.mu;

    errorP = ERROR(p_analytic, P);
    
    disp('AA');
    figure;
    contourf(C.Coll_x,C.Coll_y,P,'LineWidth',0.05) ;
    c = colorbar;
    str = {'Pressure'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Pressure','Fontsize',16)




