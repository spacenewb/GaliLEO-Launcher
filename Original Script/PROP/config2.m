close all
clear all
clc

Tin = 1.8495e+4; % N
Tfin = 5.0377e+3;
Pc = 70*10^5; % Pa
MM = 21.348;
k = 1.1841;
T0 = 2643.3; % [K]
Ru = 8314.46;
g0 = 9.81;
rho = 1.7844e+3; % [kg/m^3]
mprop = 457.1617; % [kg]
alpha = 15; % [deg]
beta = 45; % [deg]
alpha = deg2rad(alpha);
beta = deg2rad(beta);
M_c = 0.2;

[stage_2] = second_stage (Tin, Tfin, Pc, MM, k, rho, mprop)