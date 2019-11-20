clear all
clc
close all

Directory1 = 'NACA_2412_N_512.txt';
Directory2 = 'NACA_2412_N_256.txt';
Directory3 = 'NACA_2412_N_128.txt';
Directory4 = 'NACA_2412_N_64.txt';
Directory5 = 'NACA_2412_N_32.txt';
Directory6 = 'NACA_2412_N_16.txt';

AirFoil.Wing_512 = importdata(Directory1);
AirFoil.Wing_256 = importdata(Directory2);
AirFoil.Wing_128 = importdata(Directory3);
AirFoil.Wing_64 = importdata(Directory4);
AirFoil.Wing_32 = importdata(Directory5);
AirFoil.Wing_16 = importdata(Directory6);

Directory1 = 'NACA_0009_N_512.txt';
Directory2 = 'NACA_0009_N_256.txt';
Directory3 = 'NACA_0009_N_128.txt';
Directory4 = 'NACA_0009_N_64.txt';
Directory5 = 'NACA_0009_N_32.txt';
Directory6 = 'NACA_0009_N_16.txt';

AirFoil.Tail_512 = importdata(Directory1);
AirFoil.Tail_256 = importdata(Directory2);
AirFoil.Tail_128 = importdata(Directory3);
AirFoil.Tail_64 = importdata(Directory4);
AirFoil.Tail_32 = importdata(Directory5);
AirFoil.Tail_16 = importdata(Directory6);

save('AirFoil');