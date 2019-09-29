function du = diffusion_u(u,L)

N = size(u,2);
M = size(u,1);
dx = L/N;
dy = L/M;
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

for i = 2:N-1
    for j=2:M-1
        %Definim els valors de les velocitats als nodes W,E,S,N,P
        uE = u(i+1,j);
        uP = u(i,j);
        uW = u(i-1,j);
        uN = u(i,j+1);
        uS = u(i,j-1);

        du(i-1,j-1) = (uE - uP) / dx * Se + (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn + (uP-uS)/dy * Ss;
    end
end

end
