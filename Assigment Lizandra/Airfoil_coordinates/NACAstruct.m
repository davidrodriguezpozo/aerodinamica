clear all
clc
close all

Directory = 'NACA_2412_N_512.txt';
AirFoil512.Wing = importdata(Directory);
Directory = 'NACA_2412_N_256.txt';
AirFoil256.Wing = importdata(Directory);
Directory = 'NACA_2412_N_128.txt';
AirFoil128.Wing = importdata(Directory);
Directory = 'NACA_2412_N_64.txt';
AirFoil64.Wing = importdata(Directory);
Directory = 'NACA_2412_N_32.txt';
AirFoil32.Wing = importdata(Directory);
Directory = 'NACA_2412_N_16.txt';
AirFoil16.Wing = importdata(Directory);


Directory = 'NACA_0009_N_512.txt';
AirFoil512.Tail = importdata(Directory);
Directory = 'NACA_0009_N_256.txt';
AirFoil256.Tail = importdata(Directory);
Directory = 'NACA_0009_N_128.txt';
AirFoil128.Tail = importdata(Directory);
Directory = 'NACA_0009_N_64.txt';
AirFoil64.Tail = importdata(Directory);
Directory = 'NACA_0009_N_32.txt';
AirFoil32.Tail = importdata(Directory);
Directory = 'NACA_0009_N_16.txt';
AirFoil16.Tail = importdata(Directory);

clear Directory

save('AirFoil16');
save('AirFoil32');
save('AirFoil64');
save('AirFoil128');
save('AirFoil256');
save('AirFoil512');