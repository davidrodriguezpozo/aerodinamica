function du = diffusion_u(u,L)

N = size(u,2);
M = size(u,1);
dx = L/N;
dy = L/M;
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

for i = 1:N
    for j=1:M
        %Definim els valors de les velocitats als nodes W,E,S,N,P
        uE = u(i+1,j);
        uP = u(i,j);
        uW = u(i-1,j);
        uN = u(i,j+1);
        uS = u(i,j-1);

        du(i,j) = (uE - uP) / dx * Se + (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn + (uP-uS)/dy * Ss;
    end
end

end
