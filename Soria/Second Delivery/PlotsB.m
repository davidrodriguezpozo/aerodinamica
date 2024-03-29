% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function PlotsB (C,P,Vx,errorP)
    
    set(groot, 'DefaultTextInterpreter','latex');
    set(groot, 'Defaultaxesticklabelinterpreter','latex');
    
    figure; 
    loglog(1./Vx,errorP,'o-');
    title('Pressure vs. mesh size','Fontsize',18);
    xlabel('Size of the mesh h');
    ylabel('Error pressure');

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