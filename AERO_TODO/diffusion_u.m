%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%           FUNCTION DIFUSSION       %%%%%%%%%%%%
%%% This function computes the diffusive term of the    %%%
%%% NS equation.                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function du = diffusion_u(u,L)
% Inputs: N (size of the x-axis "mesh"), M (size of the y-axis "mesh"),
% dx and dy are the incrementals; and the Sx components are the Sufaces defined in the faces of the mesh. 
N = size(u,1);
M = size(u,2);
du = zeros(N,M);
dx = L/(N-2);
dy = L/(M-2);
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

for i = 2:N-1
    for j=2:M-1
        % the velocity components are definied in each face. 
        uE = u(i,j+1);
        uP = u(i,j);
        uW = u(i,j-1);
        uN = u(i+1,j);
        uS = u(i-1,j);
        % Finally, the diffusive term is computed for each point of the domain, using the values computed below. 
        du(i,j) = (uE - uP) / dx * Se - (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn - (uP-uS)/dy * Ss;
      
        
    end
end
