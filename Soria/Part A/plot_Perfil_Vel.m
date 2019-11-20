function plot_Perfil_Vel (datos,C,u,v)
    figure
    contourf(C.stagY_x,C.stagY_y,v,'LineWidth',0.05) ;
    c = colorbar;
    str = {'Y-Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Y-Velocity along the domain','Fontsize',16)
    
    figure
    contourf(C.stagX_x,C.stagX_y,u,'LineWidth',0.05) ;
    c = colorbar;
    str = {'X-Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('X-Velocity along the domain','Fontsize',16)
    
    for i = 2:datos.Nx-1
        for j = 2:datos.Ny-1
            x_v (i-1,j-1) = C.X(i);
            y_v (i-1,j-1) = C.Y(j);
            
            vel_x = (u(i,j) + u(i-1,j))/2;
            vel_y = (v(i,j) + v(i,j-1))/2;
            
            modVel (i-1,j-1) = sqrt(vel_x^2 + vel_y^2);
            
        end
    end
    
    figure
    contourf(x_v,y_v,modVel,'LineWidth',0.05) ;
    c = colorbar;
    str = {'Module of Velocity'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Module of Velocity along the domain','Fontsize',16)