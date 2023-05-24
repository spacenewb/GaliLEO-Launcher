%%% Config  %%%

%% ORBIT DATA

data.orbit.h = 400;           % orbit target altitude [km]
data.orbit.mu = 398600;       % earth gravitational constant 
data.orbit.Re = 6378;         % earth radius [km]
data.orbit.MIN_in = 45;       % Nominal inclination [°]
data.orbit.g = 9.81;          % gravity acceleration
data.orbit.rp_e = 150;        % altitude of pereigee elliptical orbit [km]
data.orbit.losses = 1200;     % dv losses for drag and gravity[m/s]

%% LAUNCHER DATA
DATA_PATH = '../DATA/';

% stage 1 
data.launcher.st1.tw_ratio = 2;                             % initial T/W ratio 
data.launcher.st1.isp = 257.3;                              % specific impulse [s]
data.launcher.st1.eps = 0.1;                                % guessed structural index                             
data.launcher.st1.C = 0.83;                                 % Diameter of the stage
data.launcher.st1.S = (data.launcher.st1.C/2)^2 * pi;       % Reference area        
data.launcher.st1.aero.full = load('full_st1.mat');         % Matrix with aerodynamic coefficient in full condition
data.launcher.st1.aero.empty = load('empty_st1.mat');       % Matric with aerodyanmic coefficient in empty condition

% stage 2
data.launcher.st2.tw_ratio = 2.0;
data.launcher.st2.isp = 290.3;
data.launcher.st2.eps = 0.1;
data.launcher.st2.aero.full = load('full_st2.mat');
data.launcher.st2.aero.empty = load('empty_st2.mat');

% stage 3
data.launcher.st3.T = [2000 2000];                          % Thrust [N]
data.launcher.st3.isp = 270;
data.launcher.st3.eps = 0.1;
data.launcher.st3.aero.full = load('full_st2.mat');
data.launcher.st3.aero.empty = load('empty_st2.mat');
data.launcher.st3.pay = 50;                                 % Payload mass [kg]
data.launcher.st3.fairing = 1;                              % Fairing length  [m]

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

data.str.E = 17e9;              % Young module of carbon fiber

%% LAUNCH SYTE 

data.site.latitude = 37.834;                                                        % Latitude of launch site 
data.site.azimuth =  asind( cosd(data.orbit.MIN_in)/cosd(data.site.latitude) );     % Launch azimuth to achieve nominal inclination
data.site.v_e = (1670/3.6) * cosd(data.site.latitude);                              % Rotating speed of Earth at laucnh site
data.site.z0 = 0;                                                                   % altitude of launch site [m]
data.site.OMEGA = 89.9;                                                             % Initial inclination of the rocket 
% wind data for the wind model
data.site.day = 150;
data.site.wind_mag = 1;
% position of launch site
data.site.lat0 = 37.834;
data.site.lon0 = -75.487;

%% PARACHUTE DATA

% Drogue 1
data.para(1).CL = 0;            % Lift coefficient
data.para(1).CD = 0.9;          % Drag coefficient 
data.para(1).S = 23.803;        % Surface of the parachute         
data.para(1).z_cut = 29902;     % Cut altitude of the parachute
% Drogue 2
data.para(2).CL = 0;
data.para(2).CD = 0.9;
data.para(2).S = 95.244;
data.para(2).z_cut = 5912;
% Main
data.para(3).CL = 0;
data.para(3).CD = 0.8;
data.para(3).S = 190.49;
data.para(3).z_cut = 0;

%% SETTINGS

% settings direct shooting  
data.optimo.N = 6; % Number of control intervals for each stage

% settings for ODE solver 
data.options = odeset('RelTol', 1e-8,'AbsTol', 1e-8,'Events',@ground);
data.options_ascent = odeset('AbsTol', 1e-8,'Events',@event_opening);
data.options_descent = odeset('AbsTol', 1e-8,'Events',@event_cut);
% settings for fsolve
data.optionsfsolve = optimoptions('fsolve','Display', 'none');
% settings for fmincon
data.optionsfmin = optimset('Algorithm','interior-point','Display','iter', 'TolX',1e-6, 'TolFun',1e-6,'TolCon',1e-6, ...
    'MaxIter', 200, 'MaxFunEvals', 100000, 'UseParallel', true);

%% Flag for different simulation
% put true or false in base of the part of the main that you want to run

data.optimize = false;
data.simulate = false;
data.control = false;
data.trajectory = true;