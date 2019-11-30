function Task22(N, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3, lt, ct1, bh, it, cdt1, cdt2, cdt3)
%% PREPROCES
% Parameters
Sw=(cw1*bw)-0.5*(bw*(cw1-cw2)); % Surface [m^2]
St=ct1*bh;
S=Sw+St;
% Coordenades dels Horsehoes (h)
angle=0;
subdivision=pi/N;
for i=1:N+1
    yhw(i)=(bw/2)*cos(angle);
    xhw(i)=0;
    yht(i)=(bh/2)*cos(angle);
    xht(i)=lt;
    yh(i)=yhw(i); 
    xh(i)=xhw(i);
    yh(N+1+i)=yht(i); 
    xh(N+1+i)=xht(i);
    angle=angle+subdivision;
end

% Coordenades dels Punts de Control (c)
c_slope=(cw1*(3/4)-(((cw1-cw2)/2) + cw2*(3/4)))/(bw/2);
for i=1:N
    ycw(i)= yhw(i)-((yhw(i)-yhw(i+1))/2);
    xcw(i)=(3/4)*cw1 - c_slope*abs(ycw(i));
    yct(i)= yht(i)-((yht(i)-yht(i+1))/2);
    xct(i)=lt+(3/4)*ct1; 
    yc(i)=ycw(i); 
    xc(i)=xcw(i);
    yc(N+i)=yct(i); 
    xc(N+i)=xct(i);
end

a=0;
for alpha=0:0.5:8
a=a+1;
% C�lculo del Twist angle wing
alphaZL=2.2*pi/180; 
tw_slope=(tw*pi/180)/(bw/2);
alpha=alpha*pi/180;
for i=1:N
    theta(i)=abs(ycw(i))*tw_slope;
    epsilon(i)=alpha+alphaZL+theta(i);
    epsilon(N+i)=alpha+(it*pi/180);
    zc(i)=0;
    zc(N+i)=-0.2*cw1;
end

for i=1:N+1
    zh(i)=0;
    zh(N+1+i)=-0.2*cw1;
end

for i=1:2*N
    nc(i,1)=sin(epsilon(i));
    nc(i,2)=0;
    nc(i,3)=cos(epsilon(i));
end


%% SOLVER
ur=[-1;0;0];
for i=1:2*N
    B(i,1)=-Uinf*nc(i);
    NC=[nc(i,1);nc(i,2);nc(i,3)];
    for j=1:2*N
        r2=[(xc(i)-xh(j+1));(yc(i)-yh(j+1)); zc(i)-zh(j+1);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infA=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[(xc(i)-xh(j));(yc(i)-yh(j));zc(i)-zh(j);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infB=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[(xc(i)-xh(j));(yc(i)-yh(j));zc(i)-zh(j);];
        r1=[(xc(i)-xh(j+1));(yc(i)-yh(j+1));zc(i)-zh(j+1);];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        AB=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        
        A(i,j)=infA+AB-infB;
    end
end

% Resoluci� de la Circulaci�
Cir=A\B;

%% POSTPROCESS
% Lift coefficient
CL=0;
delta=zeros(2*N,1);
for i=1:N
delta(i)=abs(abs(yhw(i+1))-abs(yhw(i)));
delta(N+i)=abs(abs(yht(i+1))-abs(yht(i)));
end
for i=1:N*2
CL=CL+((Cir(i)*delta(i))/(Uinf*S));
end
CL=2*CL;

% Viscous drag Coefficient
CDvw=cdw1+cdw2*CL+cdw3*CL^2;
CDvt=cdt1+cdt2*CL+cdt3*CL^2;
CDv=CDvw+CDvt;

xin=xh; zin=zh;
for i=1:N
    yin(i)=ycw(i);
    yin(N+i)=yct(i);
end

% Induced velocity
for i=1:2*N
    B(i,1)=-Uinf*nc(i,1);
    NC=[nc(i,1);nc(i,2);nc(i,3)];
    for j=1:2*N
        r2=[(xin(i)-xh(j+1));(yin(i)-(yh(j+1)));zin(i)-zh(j+1);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        VA=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec;
        
        r2=[(xin(i)-xh(j));(yin(i)-(yh(j)));zin(i)-zh(j);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        VB=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec;
        
        if i==j
            VAB=0;
        else
        r2=[(xin(i)-xh(j));(yin(i)-(yh(j)));zin(i)-zh(j);];
        r1=[(xin(i)-xh(j+1));(yc(i)-(yh(j+1)));zin(i)-zh(j+1);];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        VAB=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        end
        if mod_r1r2+r1_r2==0
        VAB=0;
        end

        Aind(i,j)=VA-VB+VAB;
    end 
    
    Vin(i)=dot(Aind(i,:),Cir);
    alpha_i(i)=Vin(i)/Uinf;
end

% Induced drag coefficient
CDi=0;
for i=1:2*N
CDi=CDi+((Cir(i)*delta(i)*alpha_i(i))/(Uinf*S));
end
CDi=-2*CDi;

% Drag coefficient
CD=CDv+CDi;

% Save results
CL_alpha(a)=CL;
CDi_alpha(a)=CDi;
CDv_alpha(a)=CDv;
CD_alpha(a)=CD;
Alpha(a)=alpha*180/pi;
end

%% RESULTS
figure(1);
plot(CL_alpha,CD_alpha,'LineWidth',1.5);
grid on;
title('Polar $C_D$-$C_L$','Fontsize',20);
ylabel('Drag coefficient $C_D$','Fontsize',16);
xlabel('Lift coefficient $C_L$','Fontsize',16);

figure(2);
plot(Alpha,CL_alpha,'r','LineWidth',1.5);
grid on;
title('$C_L$ vs $\alpha$','Fontsize',20);
xlabel('Angle of atack $\alpha$ [$^o$]','Fontsize',16);
ylabel('Lift coefficient $C_L$','Fontsize',16);

figure(3);
plot(Alpha,CD_alpha,'LineWidth',1.5); hold on;
plot(Alpha,CDi_alpha,'-.','LineWidth',1.5);
plot(Alpha,CDv_alpha,'-.','LineWidth',1.5);
grid on;
title('$C_D$ vs $\alpha$','Fontsize',20);
xlabel('Angle of atack $\alpha$ [$^o$]','Fontsize',16);
ylabel('Drag coefficient $C_D$','Fontsize',16);
legend({'Total Drag $C_D$','Induced Drag $C_{D_i}$','Viscous Drag $C_{D_v}$'}...
    ,'Interpreter','latex','Fontsize',14,'Location','Northwest')

for p = 1:length(Alpha)
    Eff(p) = CL_alpha(p)/CD_alpha(p);
end

figure(4);
plot(Alpha,Eff,'r','LineWidth',1);
grid on;
title('Efficiency $\frac{C_L}{C_D}$ vs $\alpha$','Fontsize',20);
xlabel('Angle of atack $\alpha$ [$^o$]','Fontsize',16);
ylabel('Efficiency $\frac{C_L}{C_D}$','Fontsize',16);

end

