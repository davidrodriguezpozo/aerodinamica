function [V] = semiinfinite(x2,ur,x)

R2=x-x2;
ur2=R2/(norm(R2));
ur2dotur=dot(ur,ur2);
vect=cross(ur,R2);
num=(1-ur2dotur)*vect;
denom=((norm(vect))^2)*4*pi;
V=num/denom;

end
