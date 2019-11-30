function [CL,CDv,CDi,CD] =Task21(Ground_effect, h, N, alpha, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3)
%% PREPROCES
% Parameters
S=(cw1*bw)-0.5*(bw*(cw1-cw2)); % Surface [m^2]

% Coordenades dels Horsehoes (h)
angle=0;
subdivision=pi/N;
for i=1:N+1
    yhw(i)=(bw/2)*cos(angle);
    xhw(i)=0;
    angle=angle+subdivision;
end

% Coordenades dels Punts de Control (c)
c_slope=(cw1*(3/4)-(((cw1-cw2)/2) + cw2*(3/4)))/(bw/2);
for i=1:N
    ycw(i)=(yhw(i)-((yhw(i)-yhw(i+1))/2));
    xcw(i)=(2/4)*cw1 - c_slope*abs(ycw(i));
end

% C�lculo del Twist angle wing
alphaZL=2.2*pi/180; 
tw_slope=(tw*pi/180)/(bw/2);
alpha=alpha*pi/180;
for i=1:N
    % Wing
    theta(i)=abs(ycw(i))*tw_slope;
    epsilon(i)=alpha+alphaZL+theta(i);
    nc(i,1)=sin(epsilon(i));
    nc(i,2)=0;
    nc(i,3)=cos(epsilon(i));
end

%% SOLVER
ur=[-1;0;0];
for i=1:N
    B(i,1)=-Uinf*nc(i);
    NC=[nc(i,1);nc(i,2);nc(i,3)];
    for j=1:N
        r2=[(xcw(i)-xhw(j+1));(ycw(i)-yhw(j+1));0;];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infA=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
   
        r2=[(xcw(i)-xhw(j));(ycw(i)-yhw(j));0;];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infB=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[(xcw(i)-xhw(j));(ycw(i)-yhw(j));0;];
        r1=[(xcw(i)-xhw(j+1));(ycw(i)-yhw(j+1));0;];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        AB=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        
        if Ground_effect==1
        r2=[(xcw(i)-xhw(j+1));(ycw(i)-yhw(j+1));(-2*h);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infA_ge=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[(xcw(i)-xhw(j));(ycw(i)-yhw(j));(-2*h);];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        infB_ge=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[(xcw(i)-xhw(j));(ycw(i)-yhw(j));(-2*h);];
        r1=[(xcw(i)-xhw(j+1));(ycw(i)-yhw(j+1));(-2*h);];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        AB_ge=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        
        Age = infA_ge + AB_ge - infB_ge;
        else
            Age=0;
        end
        
        A(i,j)=infA+AB-infB-Age;
    end
end

% Resoluci� de la Circulaci�
Cir=A\B;

%% POSTPROCESS

% Lift coefficient
CL=0;
for i=1:N
delta(i)=abs(abs(yhw(i+1))-abs(yhw(i)));
CL=CL+((Cir(i)*delta(i))/(Uinf*S));
end
CL=2*CL;

% Viscous drag Coefficient
CDv=cdw1+cdw2*CL+cdw3*CL^2;

% Induced velocity
yin=ycw;
for i=1:N
    B(i,1)=-Uinf*nc(i,1);
    NC=[nc(i,1);nc(i,2);nc(i,3)];
    for j=1:N
        r2=[0;(yin(i)-yhw(j+1));0;];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        VA=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec;
        
        r2=[0;(yin(i)-yhw(j));0;];
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
        r2=[0;(yin(i)-yhw(j));0;];
        r1=[0;(yin(i)-yhw(j+1));0;];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        VAB=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        end
        
        if Ground_effect==1
        r2=[0;(yin(i)-yhw(j+1));-2*h;];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        VA_ge=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        r2=[0;(yin(i)-yhw(j));-2*h;];
        mod_r2=norm(r2);
        ur2=r2/mod_r2;
        ur_ur2=dot(ur,ur2);
        urxr2=cross(ur,r2);
        mod_urxr2=norm(urxr2);
        vec=dot(urxr2,NC);
        VB_ge=(1/(4*pi))*((1-ur_ur2)/(mod_urxr2^2))*vec; 
        
        if i==j
            VAB_ge=0;
        else
        r2=[0;(yin(i)-yhw(j));-2*h;];
        r1=[0;(yin(i)-yhw(j+1));-2*h;];
        r1xr2=cross(r1,r2);
        mod_r1=norm(r1);mod_r2=norm(r2);
        mod_r1r2=mod_r1*mod_r2;
        r1_r2=dot(r1,r2);
        vec=dot(r1xr2,NC);
        VAB_ge=((mod_r1+mod_r2)/((4*pi*mod_r1r2)*(mod_r1r2+r1_r2)))*vec; 
        end
        
        V_ge = VA_ge + VAB_ge - VB_ge;
        else
            V_ge=0;
        end
        Aind(i,j)=VA-VB+VAB-V_ge;
    end 
    Vin(i)=dot(Aind(i,:),Cir);
    alpha_i(i)=Vin(i)/Uinf;
end

% Induced drag coefficient
CDi=0;
for i=1:N
CDi=CDi+((Cir(i)*delta(i)*alpha_i(i))/(Uinf*S));
end
CDi=-2*CDi;

% Drag coefficient
CD=CDv+CDi;

%% RESULTS
disp('Lift coefficient (CL): '); disp(CL);
disp('Viscous drag coefficient (CDv): '); disp(CDv);
disp('Induced drag coefficient (CDi): '); disp(CDi);
disp('Total drag coefficient (CD): '); disp(CD);

end
