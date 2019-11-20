clear all
clc
close all

A = 512;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil512.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil512.Tail = importdata(Directory);
clear Directory
save('AirFoil512');
clear all

A = 256;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil256.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil256.Tail = importdata(Directory);
clear Directory
save('AirFoil256');
clear all

A = 128;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil128.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil128.Tail = importdata(Directory);
clear Directory
save('AirFoil128');
clear all

A = 64;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil64.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil64.Tail = importdata(Directory);
clear Directory
save('AirFoil64');
clear all

A = 32;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil32.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil32.Tail = importdata(Directory);
clear Directory
save('AirFoil32');
clear all

A = 16;

Directory = sprintf('NACA_2412_N_%d.txt',A);
AirFoil16.Wing = importdata(Directory);
Directory = sprintf('NACA_0009_N_%d.txt',A);
AirFoil16.Tail = importdata(Directory);
clear Directory 
save('AirFoil16');
clear all

