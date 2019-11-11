function Plots (C, U, V, P)
   
X = C.stagX_x; 
Y = C.stagX_y;
figure
[A,h] = contourf(X,Y,U,100);
set(h,'LineColor','none')

colormap jet;
c = colorbar;
%caxis([-0.1 0.8])
str = {'X-Velocity'}; 
c.Label.String = str;
c.Label.FontSize = 16;

hold on
axis equal

xlim([ -0.1 1.1])
ylim([ -0.1 1.1])

title ('X-Velocity for Re = 100','Interpreter','latex','Fontsize',20);
xlabel('X-axis [m]','Interpreter','latex','Fontsize',16);
ylabel('Y-axis [m]','Interpreter','latex','Fontsize',16);

clear X Y
X = C.stagY_x; 
Y = C.stagY_y;
figure
[B,h2] = contourf(X,Y,V,100);
set(h2,'LineColor','none')

colormap jet;
c = colorbar;
%caxis([-0.1 0.8])
str = {'Y-Velocity'}; 
c.Label.String = str;
c.Label.FontSize = 16;

hold on
axis equal

xlim([ -0.1 1.1])
ylim([ -0.1 1.1])

title ('Y-Velocity for Re = 100','Interpreter','latex','Fontsize',20);
xlabel('X-axis [m]','Interpreter','latex','Fontsize',16);
ylabel('Y-axis [m]','Interpreter','latex','Fontsize',16);

clear X Y

X = C.Coll_x; 
Y = C.Coll_y;
figure
[D,h3] = contourf(Y,X,P,100);
set(h3,'LineColor','none')

colormap jet;
c = colorbar;
%caxis([-0.1 0.8])
str = {'Pressure'}; 
c.Label.String = str;
c.Label.FontSize = 16;

hold on
axis equal

xlim([ -0.1 1.1])
ylim([ -0.1 1.1])

title ('Pressure for Re = 100','Interpreter','latex','Fontsize',20);
xlabel('X-axis [m]','Interpreter','latex','Fontsize',16);
ylabel('Y-axis [m]','Interpreter','latex','Fontsize',16);