function [P_new] = pressureC(datos, u, v, u_p, v_p, delta_t,P)
%Solve the system Ap = b
delta = datos.L/datos.Vx;
P_new  = zeros(datos.Nx);
a_e = zeros(datos.Nx);
a_w = zeros(datos.Nx);
a_s = zeros(datos.Nx);
a_n = zeros(datos.Nx);
a_p = zeros(datos.Nx);
b_p = zeros(datos.Nx);
    
for i = 2:datos.Vx+1
    for j = 2:datos.Vx+1
    a_e(i,j) = 1;
    a_w(i,j) = 1;
    a_s(i,j) = 1;
    a_n(i,j) = 1;
    a_p(i,j) = 4;
    b_p(i,j) = datos.rho/delta_t*(u_p(i,j)*delta-u_p(i,j-1)*delta+v_p(i,j)*delta-v_p(i-1,j)*delta);
    end
end
P_prev = P;
error = 10;

while error > 0.1
    for i = 2:datos.Vx+1
        for j = 2:datos.Vx+1
            P(i,j) = (a_e(i,j)*P(i,j+1)+a_w(i,j)*P(i,j-1)+a_s(i,j)*P(i-1,j)+a_n(i,j)*P(i+1,j)+b_p(i,j))/a_p(i,j);
        end
    end
    
P = haloupdate(P);

error_i = 0;

for i=2:datos.Nx-1
    for j = 2:datos.Nx-1
        if abs(P(i,j) - P_prev(i,j)) > error_i
            error_i = abs(P(i,j) - P_prev(i,j));
        end
    end
end
error = error_i;
P_prev = P;
end
