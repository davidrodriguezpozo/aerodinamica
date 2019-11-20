function [cp] = Coefficients_Pressure(d) 

disp ('Calculo de coeficientes presi�n');

ae = ones(d.Nrow, d.Ncol);
aw = ones(d.Nrow, d.Ncol);
an = ones(d.Nrow, d.Ncol);
as = ones(d.Nrow, d.Ncol);
bp = zeros(d.Nrow, d.Ncol);

% FRONTERAS

aw(:,1) = zeros(d.Nrow,1);
an(:,1) = zeros(d.Nrow,1);
as(:,1) = zeros(d.Nrow,1);

ae(:,d.Ncol) = zeros(d.Nrow,1);
an(:,d.Ncol) = zeros(d.Nrow,1);
as(:,d.Ncol) = zeros(d.Nrow,1);

for i = 1:d.Nrow-1
    ae(i,1) = 2; aw(i,2) = 2;
    aw(i,d.Ncol) = 2; ae(i,d.Ncol-1) = 2;
end

ae(1,:) = zeros(1,d.Nrow);
aw(1,:) = zeros(1,d.Nrow);
an(1,:) = zeros(1,d.Nrow);

as(d.Nrow,:) = zeros(1,d.Nrow);
aw(d.Nrow,:) = zeros(1,d.Nrow);
ae(d.Nrow,:) = zeros(1,d.Nrow);

for j = 1:d.Ncol-1
    as(1,j) = 2; an(2,j) = 2;
    an(d.Nrow,j) = 2; as(d.Nrow-1,j)=2;
end

ap = ae + aw + an + as;

cp.ae = ae;
cp.aw = aw;
cp.an = an;
cp.as = as;
cp.ap = ap;
cp.bp = bp;

disp ('Fin Calculo de coeficientes presi�n');