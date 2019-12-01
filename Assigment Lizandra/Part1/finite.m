function V = finite(x1,x2,x)
r1=x-x1;
r2=x-x2;
R1=norm(r1);
R2=norm(r2);
r1dotr2=dot(r1,r2);
num=R1+R2;
denom=R1*R2*(R1*R2+r1dotr2);
V=cross(r1,r2);
V=V*(num/denom);
V=V/(4*pi);
end

