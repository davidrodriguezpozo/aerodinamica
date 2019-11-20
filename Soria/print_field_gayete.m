function []=print_field(u,name)
N=size(u,2);
M=size(u,1);

%%for i=1:lenght(u(:,1))
for j=M:-1:1
    fprintf('j=%d ; \t',j);
    for i=1:N
       fprintf('%.2f . ',u(i,j));
    end
    fprintf('\n');
end



end
