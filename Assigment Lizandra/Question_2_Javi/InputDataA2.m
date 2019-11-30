%% INPUT DATA MAIN A2
%% General data
alpha=4;   % Angle d'atac [·]
lt=3.2;     % Separació Ala/Cua [m]
h=1;        % Altura [m]
N=64;       % Segments de ala
Uinf=40;
%% WING: NACA 2412 (Trapezoidal)
NACAtail='NACA_2412_N_16.txt';
cw1=1.2;    % Corda arrel-ala [m]
cw2=0.8;    % Corda punta-ala [m]
bw=6;       % Emvergadura [m]
tw=-3;      % Twist [·]
% Viscous drag coefficient: Cdvw=cdvw1+cdvw2*Cl+cdvw3*Cl^2
cdw1=0.0063;
cdw2=-0.0033;
cdw3=0.0067;

%% TAIL: NACA 0009 (Rectangular)
NACAtail='NACA_0009_N_16.txt';
ct1=0.6;    % Corda cua [m]
bh=3;       % Envergadura cua [m]
it=-3;      % Angle de incidència [·]
cdt1=0.0055;
cdt2=0;
cdt3=0.0045;

%% Importar data dels perfils
wing=fopen ('NACA_2412_N_16.txt','r');
Dataw=textscan(wing,'%f %f %f');
fclose (wing);
dataw2=Dataw{2};
dataw3=Dataw{3};

tail=fopen (NACAtail,'r');
Datat=textscan(tail,'%f %f %f');
fclose (tail);
datat2=Datat{2};
datat3=Datat{3};

