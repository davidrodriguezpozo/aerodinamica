function du = diffusion_u(u,L)

N = size(u,2);
M = size(u,1);
du = zeros(N,M);
dx = L/(N-2);
dy = L/(M-2);
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn; 

for i = 2:N-1
    for j=2:M-1
        %Definim els valors de les velocitats als nodes W,E,S,N,P
        uE = u(i,j+1);
        uP = u(i,j);
        uW = u(i,j-1);
        uN = u(i+1,j);
        uS = u(i-1,j);

        %du(i,j) = (uE - uP) / dx * Se + (uP-uW)/ dx * Sw + (uN-uP)/dy * Sn + (uP-uS)/dy * Ss;
        
        dx = 0.2;
      
        du(i,j)=(u(i,j+1)-u(i,j))/dx*dx-(u(i,j)-u(i,j-1))/dx*dx...
          + (u(i+1,j)-u(i,j))/dx*dx-(u(i,j)-u(i-1,j))/dx*dx;
        a = 1;
    end
end