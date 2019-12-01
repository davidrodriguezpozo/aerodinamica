function [sa ca li] = compute_angles (x1,x2,z1,z2)
li = sqrt((x1-x2)^2+(z1-z2)^2);
sa = (z1-z2) / (li); 
ca = (x2-x1) / (li);
end