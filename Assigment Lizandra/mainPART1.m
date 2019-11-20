clear all
close all
clc

ind = 16;

AirFoil_name = sprintf ('Airfoil%d',ind);

Airfoil = load(AirFoil_name);
disp('Airfoil coordinates');
