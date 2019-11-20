function cu =convection(u,v,L)
%Cal diferenciar entre velocitat horitzontal i velocitat vertical. 
M=size(u,1);
N=size(v,2);
dx = L/N;
dy = L/M;
Sw = dy;
Sn = dx;
Se = Sw;
Ss = Sn;


for i=1:N
    for j=1:M
         %Definim els valors de les velocitats als nodes W,E,S,N,P
        uP =  u(i,j);
        uE = (u(i+1,j)+uP)/2;
        uW = (u(i-1,j)+uP)/2;
        uN = (u(i,j+1)+uP)/2;
        uS = (u(i,j-1)+uP)/2;
        
        %Definim els valors del fluxe 
        
        FS = dx*(v(i,j-1)+v(i+1,j-1))/2;
        FN = dx*(v(i,j)+v(i+1,j))/2;
        FE = dy*(u(i+1,j)+uP)/2;
        FW = dy*(u(i-1,j)+uP)/2;
        
        %La convecció és:
        
        
        cu(i,j) = uE*FE-uW*FW+uN*FN-uS*FS;
    end
end

end
