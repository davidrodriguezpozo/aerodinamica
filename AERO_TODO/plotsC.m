function plotsC (u, v, P, P_an, P_num, u_an, u_num, v_an, v_num, acu_time, C, ii, jj)

% 
figure(1); 
plot(acu_time,P_an); hold on;
plot(acu_time,P_num);
legend('Analytic','Numerical');
posx = num2str(C.stagX_x(ii,jj));
posy = num2str(C.stagX_y(ii,jj));
title(['Pressure at point x =',posx,' and y = ',posy]);
xlabel('time [s]');
ylabel('Pressure');

figure(2); 
plot(acu_time,u_an); hold on;
plot(acu_time,u_num);
legend('Analytic','Numerical');
posx = num2str(C.stagX_x(ii,jj));
posy = num2str(C.stagX_y(ii,jj));
title(['Horizontal Velocity at point x =',posx,' and y = ',posy]);
xlabel('time [s]');
ylabel('X-Velocity ');

figure(3); 
plot(acu_time,v_an); hold on;
plot(acu_time,v_num);
legend('Analytic','Numerical');
posx = num2str(C.stagY_x(ii,jj));
posy = num2str(C.stagY_y(ii,jj));
title(['Vertical Velocity at point x =',posx,' and y = ',posy]);
xlabel('time [s]');
ylabel('Y-Velocity ');

figure(4);
    contourf(C.Coll_x,C.Coll_y,P,'LineWidth',0.05) ;
    c = colorbar;
    str = {'Pressure'}; 
    c.Label.String = str;
    c.Label.FontSize = 14;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Pressure','Fontsize',16)