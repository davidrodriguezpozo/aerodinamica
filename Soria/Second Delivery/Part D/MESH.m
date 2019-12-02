function C = MESH (d)

fprintf('-( 1 )- Calculando mallas...\n');

% Preallocation of variables

C.X = zeros (d.Nrow,d.Ncol);
C.Y = zeros (d.Nrow,d.Ncol);

C.dx = zeros (d.Nrow,1);
C.dy = zeros (d.Ncol,1);

C.stagX_x = zeros (d.Nrow_sx,d.Ncol_sx);
C.stagX_y = zeros (d.Nrow_sx,d.Ncol_sx);

C.stagY_x = zeros (d.Nrow_sy,d.Ncol_sy);
C.stagY_y = zeros (d.Nrow_sy,d.Ncol_sy);

[ C.X, C.Y, C.dx, C.dy, C.Coll_x, C.Coll_y ] = COLLOCATED_MESH ( d );

[ C.stagX_x, C.stagX_y ] = STAGG_MESH_X ( d, C );

[ C.stagY_x, C.stagY_y ] = STAGG_MESH_Y ( d, C );

fprintf('---> Mallas calculada\n')


function [ X, Y, dx, dy, Coll_X, Coll_Y ] = COLLOCATED_MESH ( d )

% Delta X

dx = 1/d.Nx;
dy = 1/d.Nx;

C.dx = dx;
C.dy = dy;

% COLLOCATED MESHES

X (1) = 0; X(2) = dx/2;
for i = 3:d.Ncol-1
    X(i) = X(i-1) + dx;
end
X(d.Ncol) = X(d.Ncol-1) + dx/2;

Y (1) = 0; Y(2) = dy/2;
for j = 3:d.Nrow-1
    Y(j) = Y(j-1) + dy;
end
Y(d.Nrow) = Y(d.Nrow-1) + dy/2;

for i=1:d.Nrow
    for j=1:d.Ncol
        Coll_X(i,j) = abs(d.L-X(i));
        Coll_Y(i,j) = Y(j);%abs(d.L-Y(j)); % Y(j);
    end
end

fprintf(' Collocated mesh calculada\n');


function [ x, y ] = STAGG_MESH_X ( d, C )

vector_x = zeros( d.Nx ,1 );
vector_x(1) = 0; % La malla comen�a a la paret
for i = 2:d.Nx
    vector_x(i) = vector_x(i-1) + C.dx;
end
vector_x(d.Nx+1) = d.L;

for i = 1:d.Nrow_sx
    for j = 1:d.Ncol_sx
        x(i,j) = vector_x(j);
        y(i,j) = abs(d.L-C.Y(i+1));
    end
end

fprintf(' Staggered Mesh X calculada\n');


function [ X, Y ] = STAGG_MESH_Y ( d, C )

vector_y = zeros( d.Ncol_sy ,1 );
vector_y(1) = 0; % La malla comen�a a la paret

for i = 2:d.Ny
    vector_y(i) = vector_y(i-1) + C.dy;
end
vector_y(d.Ny+1) = d.L;

for i = 1:d.Nrow_sy
    for j = 1:d.Ncol_sy
        X(i,j) = C.X(j+1);
        %Y(i,j) = vector_y(i);
        Y(i,j)= vector_y(d.Nrow_sy+1-i);
    end
end

fprintf(' Staggered Mesh Y calculada\n');