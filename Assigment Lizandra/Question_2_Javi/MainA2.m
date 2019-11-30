%% ASSIGNMENT 2: AERODYNAMICS
% Asignatura de "AERODIN�MICA, MEC�NICA DE VOL I ORBITAL"
% Javier Jos� Monerris Valenti
% Marcel Mar�n de Yzaguirre
% Arnau Verdaguer Pons
% Pau Nadal Vila

% 06/12/2019

clc; clear all; close all;

%% INPUT DATA
run('InputDataA2');

%% QUESTION 1

%% QUESTION 2

set(groot, 'DefaultTextInterpreter','latex');
set(groot, 'Defaultaxesticklabelinterpreter','latex');

N = 128;
% Task 2.1:
disp('GROUND EFFECT RESULTS:')
Ground_effect=1;
Task21(Ground_effect, h, N, alpha, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3);

disp('NO GROUND EFFECT RESULTS:')
Ground_effect=0;
Task21(Ground_effect, h, N, alpha, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3);

% Task 2.2:
Task22(N, Uinf, cw1, cw2, bw, tw, cdw1, cdw2, cdw3, lt, ct1, bh, it, cdt1, cdt2, cdt3);
