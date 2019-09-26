%function PrintFacil(u, 'name')

u = [11 12 13
    21 22 23
    31 32 33];

uNew = zeros(size(u));

for i= 1:size(u, 1)
    for j=i:size(u,2)
        ind = N-i-1;
        uNew(i,j) = u(ind, j);
    end
end
        
%fprintf('hi %.2f \n', uNew)



%return 'name' = 
