% Code developed by:
% - Sergi Martinez Castellarnau
% - Carlos Perez Ricardo
% - David Rodriguez Pozo
% - Paula Sorolla Bayod

function B=haloupdate(B)

N=size(B,1)-2;
M=size(B,2)-2;
B(1,:)=B(N+1,:);
B(N+2,:)=B(2,:);
B(:,1)=B(:,M+1);
B(:,M+2)=B(:,2);

end