%%% Config  %%%

%% ORBIT DATA

data.orbit.h = 400;           % km
data.orbit.mu = 398600;       % earth gravitational constant
data.orbit.Re = 6378;         % earth radius km
data.orbit.MIN_in = 45;       % 
data.orbit.g = 9.81;
data.orbit.mu = 398600;
data.orbit.rp_e = 150;
data.orbit.losses = 1200; %[m/s]



%% LAUNCHER DATA
DATA_PATH = '../DATA/';

% stage 1 
data.launcher.st1.tw_ratio = 2; % 2.5;
data.launcher.st1.isp = 257.3;
data.launcher.st1.eps = 0.1574;
data.launcher.st1.C = 0.83;
data.launcher.st1.S = (data.launcher.st1.C/2)^2 * pi;
data.launcher.st1.aero.full = load('full_st1.mat');
data.launcher.st1.aero.empty = load('empty_st1.mat');

% stage 2
data.launcher.st2.tw_ratio = 2.0;
data.launcher.st2.isp = 290.3;
data.launcher.st2.eps = 0.1809;
data.launcher.st2.aero.full = load('full_st2.mat');
data.launcher.st2.aero.empty = load('empty_st2.mat');

% stage 3
data.launcher.st3.T = [2000 2000];
data.launcher.st3.isp = 270;
data.launcher.st3.eps = 0.2210;
data.launcher.st3.aero.full = load('full_st2.mat');
data.launcher.st3.aero.empty = load('empty_st2.mat');
data.launcher.st3.pay = 50; %[kg]
data.launcher.st3.fairing = 1; %[m]

%% PROPELLANT DATA

data.prop.st1.Pc = 70*10^5; % Pa
data.prop.st1.MM = 21.348;
data.prop.st1.k = 1.1841;
data.prop.st1.T0 = 2643.3; % [K]
data.prop.st1.Ru = 8314.46;
data.prop.st1.rho = 1.7844e+3; % [kg/m^3]
data.prop.st1.alpha = deg2rad(15); % [deg]
data.prop.st1.beta = deg2rad(45); % [deg]
data.prop.st1.Pe = 101325;
data.prop.st1.M_c = 0.2;

data.prop.st2.Pc = 70*10^5; % Pa
data.prop.st2.MM = 21.348;
data.prop.st2.k = 1.1841;
data.prop.st2.T0 = 2643.3; % [K]
data.prop.st2.Ru = 8314.46;
data.prop.st2.rho = 1.7844e+3; % [kg/m^3]
data.prop.st2.alpha = deg2rad(15); % [deg]
data.prop.st2.beta = deg2rad(45); % [deg]
data.prop.st2.M_c = 0.2;
data.prop.st2.eps = 40;
data.prop.st2.Pamb = 0;

data.prop.st3.Pc = 10*10^5; % Pa
data.prop.st3.MM = 22.256;
data.prop.st3.Cp = 2.2125;
data.prop.st3.Ru = 8314.46;
data.prop.st3.Cv = data.prop.st3.Cp - data.prop.st3.Ru/(data.prop.st3.MM*1000);
data.prop.st3.k = data.prop.st3.Cp/data.prop.st3.Cv;
data.prop.st3.T0 = 2131.11; % [K]
data.prop.st3.rho = 1357; % [kg/m^3]
data.prop.st3.alpha = deg2rad(15); % [deg]
data.prop.st3.beta = deg2rad(45); % [deg]
data.prop.st3.eps = 200;
data.prop.st3.Pamb = 0;
data.prop.st3.M_c = 0.2;

%% STRUCTURAL DATA

data.str.FOS = 2;
data.str.sigma = 110e6;
data.str.sigma_hoop = 110e6; 
data.str.E = 17e9;
data.str.rho_shell = 1840;
data.str.rho_ins = 1180;


%% LAUNCH SYTE 

data.site.latitude = 37.834;
data.site.azimuth =  asind( cosd(data.orbit.MIN_in)/cosd(data.site.latitude) );
data.site.v_e = (1670/3.6) * cosd(data.site.latitude);
data.site.z0 = 0;
data.site.OMEGA = 89.9;
data.site.day = 150;
data.site.wind_mag = 1;
data.site.lat0 = 37.834;
data.site.lon0 = -75.487;

%% PARACHUTE DATA

data.para(1).CL = 0;
data.para(1).CD = 0.9;
data.para(1).S = 23.803;
data.para(1).z_cut = 29902;
data.para(2).CL = 0;
data.para(2).CD = 0.9;
data.para(2).S = 95.244;
data.para(2).z_cut = 5912;
data.para(3).CL = 0;
data.para(3).CD = 0.8;
data.para(3).S = 190.49;
data.para(3).z_cut = 0;

%% SETTINGS

% optimo settings 
data.optimo.N = 6; % ottimizazione ok con questi parametro 6 


% option for matlab function 
data.options = odeset('RelTol', 1e-8,'AbsTol', 1e-8,'Events',@ground);
data.options_ascent = odeset('AbsTol', 1e-8,'Events',@event_opening);
data.options_descent = odeset('AbsTol', 1e-8,'Events',@event_cut);
data.optionsfsolve = optimoptions('fsolve','Display', 'none');
data.optionsfmin = optimset('Algorithm','interior-point','Display','iter', 'TolX',1e-6, 'TolFun',1e-6,'TolCon',1e-6, ...
    'MaxIter', 100, 'MaxFunEvals', 100000, 'UseParallel', true);

% checks for plot 
data.plot_opt = false;
data.plot_traj = true;
data.optimize = false;
data.simulate = false;
data.simulate_final = false;
data.recovery = true;
