function delta_T = TimeStep (datos,u,v)
        %Time for convective term
    Tc_u = min(datos.L/(datos.Vx*max(max(u))));
    Tc_v = min(datos.H/(datos.Vy*max(max(v))));
    Tc = min (Tc_u,Tc_v);
    %Time for diffusive term
    Td = 0.5*(datos.L/datos.Vx)*(datos.H/datos.Vy) / datos.mu;
    
    %Step time:
    delta_T =0.2 * min(Tc,Td);
    delta_T =0.02 ;