function [sol] = traj_opt(launcher,data)

site = data.site;
orbit = data.orbit;

% set guess vector for optimization variable
N = data.optimo.N;
tf_guess = launcher.st3.tb - 5;         % time burning of third stage as guess 
u_guess =  deg2rad(0).*ones(1,(3*N));   % Control guess TVC angle (set = 0°)

% target definition 

% evaluation of orbit parameter
rp = orbit.rp_e + orbit.Re;     % radius of perigee transfer orbit
ra = orbit.h + orbit.Re;        % radius of apogee transfer orbit 
e = (ra-rp)/(ra+rp);            % eccentricity 
a = (ra+rp)/2;                  % semi-major axis 
p = a*(1-e^2);                  % semi latus rectum
vp = sqrt(orbit.mu/p)*(1+e);    % velocity at perigee 

target.z = orbit.rp_e * 1e3;    % target altitude 
target.vx = vp * 1000;          % target velocity 
target.vz = 0;                  % target vertical velocity (0 because at perigee)

% set initial condition 
x0 = [tf_guess, u_guess];
% set boundaries
lb = [1,-deg2rad(6.5)*ones(1,(3*N))];
ub = [launcher.st3.tb, deg2rad(6.5)*ones(1,(3*N))];

% cost function 
fun = @(x) dv_eval(x,launcher,orbit,site,N,data.options);

% non-linear constraint
constraint = @(x) nonlcon(x,launcher, site, N,target, orbit, data.options);

% optimization algorithm 
[sol, DV] = fmincon(fun, x0, [], [], [], [], lb, ub, constraint, data.optionsfmin);
