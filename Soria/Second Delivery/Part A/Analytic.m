%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION ANALYTIC       %%%%%%%%%%%%%
%%% This function allows the programer to compute an    %%%
%%% analytical solution in order to validate the code   %%%
%%% afterwards. In this part, a periodic and divergen-  %%%
%%% ce-free velocity fiel is arranged.                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [conv diff] = Analytic (datos, C, malla)
% The Outputs: the convective and diffusive terms. 
% The Inputs: the mesh (C), the specific data for solving the problem and the "malla" field,
% which indicates the name of the mesh

    for i = 2:datos.Nx-1
        for j=2:datos.Ny-1
        % Firstly Delta X, Delta Y and the Surface are obtained
        dx = C.dx(i);
        dy = C.dy(j);
        Sup = dx*dy;
        
        % The program chooses if the mesh is the x-axis staggered mesh or the y-axis staggered mesh.
        % Thus, one of the staggered meshes is used.
        if strcmp(malla,'x') == 1
            x = C.stagX_x(i,j);
            y = C.stagX_y(i,j);
        elseif strcmp(malla,'y') == 1
            x = C.stagY_x(i,j);
            y = C.stagY_y(i,j);
        end
        
        % The velocity field is defined by two known velocities
        u = datos.F*cos(2*pi*x)*sin(2*pi*y);
        v = -datos.F*cos(2*pi*y)*sin(2*pi*x);
        
        % Derivatives from that field are defined (regarding both variables)
        du_dx = - datos.F*2*pi*sin(2*pi*x)*sin(2*pi*y);
        du_dy = datos.F*2*pi*cos(2*pi*x)*cos(2*pi*y);
        dv_dx = -datos.F*2*pi*cos(2*pi*y)*cos(2*pi*x);
        dv_dy = datos.F*2*pi*cos(2*pi*y)*sin(2*pi*x);
        
        % Second derivatives are computed
        d2u_dx2 = - datos.F^2*4*pi^2*sin(2*pi*y)*cos(2*pi*x);
        d2v_dy2 = datos.F^2*4*pi^2*sin(2*pi*x)*cos(2*pi*y);
        
        % Like steps before, the code asks if it is the x-axis mesh or the y-axis. 
        % Then and finally, the convective and diffusive terms are computed for each mesh. 
        if strcmp(malla,'x') == 1
            
            conv(i-1,j-1) = Sup*(u*du_dx + u*dv_dx); %u * du/dx + u * dv/dx
            diff(i-1,j-1) = Sup*(d2u_dx2); 
            
        elseif strcmp(malla,'y') == 1
            
            conv(i-1,j-1) = Sup*(v*du_dy + v*dv_dy); %v * du/dy + v * dv/dy
            diff(i-1,j-1) = Sup*(d2v_dy2); 
        
        end 
    end
    end
