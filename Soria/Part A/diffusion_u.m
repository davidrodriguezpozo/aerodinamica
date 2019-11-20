  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         FUNCTION DIFFUSION        %%%%%%%%%%%%%
%%% This function computes the diffusive term, using the%%%
%%% velocity field. Besides, it is used the notation of %%%
%%% relative coordenates from the center of the VC.     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function du = diffusion_u(u,L)
% The inputs: velocity field (u) and initial conditions (L).
% The outputs: the diffusive term of the velocity field (u or v, depending on the inputs assigned).

% N and M are the sizes of the velocity field (or the meshes). 
N = size(u,2);
M = size(u,1);

% dx and dy are the alpha X and alpha Y, the incrementals (it is assumed uniform mesh).
dx = L/N;
dy = L/M;

% The surface and the diferentials are equals because of the definition of the supposed unitary depth of the mesh. 
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

% The bucle starts in i=2 because it is computed for the i-1 point in the mesh. 
for i = 2:N-1
    for j=2:M-1
        % The velocity values at W,E,S,N,P nodes are computed. 
        uE = u(i+1,j);
        uP = u(i,j);
        uW = u(i-1,j);
        uN = u(i,j+1);
        uS = u(i,j-1);
        % Finally, the diffusive term is computed in the i-1 and j-i position, using the values previously computed. 
        du(i-1,j-1) = (uE - uP) / dx * Se + (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn + (uP-uS)/dy * Ss;
    end
end

end
