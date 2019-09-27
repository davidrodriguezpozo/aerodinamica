function []=halo_update(u)
N=size(u,1)-2;
% Columna 1
fprintf('%s \n','Columna 1');
for i=2:N+1
    u(i,1)=u(i,N+1);
    fprintf('%.2f . ',u(i,1));
end
fprintf('\n');

fprintf('%s \n','Fila 1');
for i=2:N+1
    u(1,i)=u(N+1,i);
    fprintf('%.2f . ',u(1,i));
end
fprintf('\n');

fprintf('%s \n','Columna final ');
for i=2:N+1
    u(i,N+2)=u(i,2);
    fprintf('%.2f . ',u(i,N+2))
end

fprintf('\n');
fprintf('%s \n','Fila final ');
for i=2:N+1
    u(N+2,i)=u(2,i);
    fprintf('%.2f . ',u(N+2,i))
end
fprintf('\n');


end

