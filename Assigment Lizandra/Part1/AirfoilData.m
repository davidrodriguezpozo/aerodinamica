% Airfoil data for Assignament 2.
%   Barrachina, Victor
%   El Jarari, Younes
%   Royo, Enric

% Physic parameters
rho = 1.225; %[kg/m3]
Uinf = 1;

% Main wing
c_w1 = 1.2; %[m]
c_w2 = 0.8; %[m]
b_w = 6; %[m]
C_D0_w = 0.0063; %[-]
C_D1_w = 0.0033; %[-]
C_D2_w = 0.0067; %[-]
theta_w = -3*pi/180; %[rad]
i_w = 0; %[rad]
NACA_w = [2 4 12];

% Tail wing
c_t = 0.6; %[m]
b_t = 3; %[m]
C_D0_t = 0.0055; %[-]
C_D1_t = 0.0045; %[-]
theta_t = 0; %[rad]
i_t = -3*pi/180; %[rad]
NACA_t = [0 0 09];

% Wing-tail
lt = 3.2; %[m]