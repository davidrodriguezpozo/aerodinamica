function [u_p, v_p, R_u, R_v] = predicted(u,v,delta_t,datos, conv_u, diff_u, conv_v, diff_v, R_uant, R_vant,R_u,R_v)
 
Vx = datos.Vx;
Vy = datos.Vy;
Vol = (datos.L/datos.Vx)^2;
u_p = zeros(Vx+2,Vy+2);
v_p = zeros(Vx+2,Vy+2);
for i=2:Vx+1
    for j=2:Vy+1

    R_u(i,j) = (-conv_u(i,j) + datos.mu*diff_u(i,j))/Vol;
    R_v(i,j) = (-conv_v(i,j) + datos.mu*diff_v(i,j))/Vol;
    u_p(i,j) = u(i,j)+delta_t*(3/2*R_u(i,j)+1/2*R_uant(i,j));
    v_p(i,j) = v(i,j)+delta_t*(3/2*R_v(i,j)+1/2*R_vant(i,j));
    
    end
end


end