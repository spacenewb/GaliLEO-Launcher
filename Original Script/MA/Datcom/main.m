%{
CONFIG - This script sets up all the parameters for missile Datcom
All the parameters are stored in the "datcom" structure.

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 18/10/2019

%}
clear 
close all
clc

%% States
% State values in which the aerodynamic coefficients will be computed
datcom.Mach = [0.1 0.5 0.8 1 1.2 2 4 5 6 8 10 12 14 16 18 20];
datcom.Alpha = [-20 -10 -7.5 -5 -2.5 -1.5 -1 -0.5 -0.1 0 0.1 0.5 1 1.5 2.5 5 7.5 10 20 ];
datcom.Beta = 0;
datcom.Alt = [0 500 1000 2000 3000 5000 8000 10000 12000 14000 16000 20000 22000 24000 26000 28000 30000 35000 40000 50000] ;

%% Design Parameters
datcom.Chord1 = 0.10; 
datcom.Chord2 = 0.8; 
datcom.Height = 0.10;                            
datcom.shape = 'rect';

%% Fixed Parameters
vars.xcg = [2,1];                          % [m] CG position [full, empty]
datcom.D = 0.8;                                 % [m] rocket diameter
datcom.Lnose = 1.5;                             % [m] nose length
datcom.Lcenter = 12-1.5;                        % [m] Lcenter : Centerbody length
datcom.Npanel = 4;                              % [m] number of fins
datcom.Phif = [0 90 180 270];                      % [deg] Angle of each panel
datcom.Ler = 0.003;                             % [deg] Leading edge radius
datcom.d = 0;                                   % [m] rocket tip-fin distance
datcom.zup_raw = 0.0015;                        % [m] fin semi-thickness 
datcom.Lmaxu_raw = 0.006;                       % [m] Fraction of chord from leading edge to max thickness
%% Ogive parameters
datcom.OgType = 'KARMAN';
datcom.NosePower = 1/2;

%% Protuberance parameters
datcom.xprot = 1; % axial position 
datcom.nloc = 3; % number of brakes
datcom.lprot = 0.005; % brakes thickness
datcom.wprot = 0.116; % brakes width
vars.hprot = 0; % brakes length, first entry must be always 0!

%% Run 
autoMatricesProtub(datcom, vars);
