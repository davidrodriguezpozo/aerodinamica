function V = horseshoe(x1,x2,ur,x)

v1=semiinfinite(x1,ur,x);
v12=finite(x1,x2,x);
v2=semiinfinite(x2,ur,x);
V=v1+v12-v2;

end
