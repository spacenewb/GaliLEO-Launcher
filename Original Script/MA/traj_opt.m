function [sol] = traj_opt(launcher,data)

site = data.site;
orbit = data.orbit;

% %% Adimensionalize 
% 
% param.z_ref = 400e3;
% param.t_ref = launcher.st1.tb + launcher.st2.tb + launcher.st3.tb;
% param.v_ref = param.z_ref/param.t_ref;
% param.a_ref = param.z_ref/param.t_ref.^2;
% param.m_ref = launcher.st1.m0;
% param.N_ref = param.m_ref * param.a_ref;
% 
% 
% % %% DIRECT SHOOTING 
% % 
% % % adimension parameter 
% site.z0 = site.z0/param.z_ref;
% orbit.Re = orbit.Re/param.z_ref;
% orbit.g = orbit.g/param.a_ref;
% launcher.st1.tb = launcher.st1.tb/param.t_ref;
% launcher.st2.tb = launcher.st2.tb/param.t_ref;
% launcher.st3.tb = launcher.st3.tb/param.t_ref;
% launcher.st1.mp = launcher.st1.mp/param.m_ref;
% launcher.st2.mp = launcher.st2.mp/param.m_ref;
% launcher.st3.mp = launcher.st3.mp/param.m_ref;
% launcher.st1.m0 = launcher.st1.m0/param.m_ref;
% launcher.st2.m0 = launcher.st2.m0/param.m_ref;
% launcher.st3.m0 = launcher.st3.m0/param.m_ref;
% launcher.st1.T = launcher.st1.T/param.N_ref;
% launcher.st2.T = launcher.st2.T/param.N_ref;
% launcher.st3.T = launcher.st3.T/param.N_ref;

% set guess vector for optimization variable
N = data.optimo.N;
tf_guess = launcher.st3.tb - 5; % ottimizazione ok con questo parametro 
u_guess =  deg2rad(0).*ones(1,(3*N));

% target definition 
rp = orbit.rp_e + orbit.Re;
ra = orbit.h + orbit.Re;
e = (ra-rp)/(ra+rp);
a = (ra+rp)/2;
p = a*(1-e^2);
vp = sqrt(orbit.mu/p)*(1+e);

target.z = orbit.rp_e * 1e3;
target.vx = vp * 1000;
target.vz = 0;

% set boundaries and intial condition 
x0 = [tf_guess, u_guess];
lb = [1,-deg2rad(6.5)*ones(1,(3*N))];
ub = [launcher.st3.tb, deg2rad(6.5)*ones(1,(3*N))];

fun = @(x) dv_eval(x,launcher,orbit,site,N,data.options);
constraint = @(x) nonlcon(x,launcher, site, N,target, orbit, data.options);

[sol, DV] = fmincon(fun, x0, [], [], [], [], lb, ub, constraint, data.optionsfmin);


% sol(1) = sol(1)*param.t_ref





% %% DIRECT MULTIPLE  SHOOTING 
% 
% N = data.optimo.N;
% tf_guess = 1000;
% u_guess =  0.*ones(N+1,1);
% 
% x_initial = zeros(1,8); % number of states x , y , vx ecc
% x_in_guess = zeros(N+1,8);
% h = tf_guess/N;
% 
% for i = 1:N
%     
%     [~,Y] = ode45(@rocket_dyn,[(i-1)*h h*i],x_initial,options);
%     x_initial = Y(end,:);
%     x_in_guess = x_sol(end,:);
%     
% end
% 
% x_guess_r = reshape(x_guess,[],1);
%   
% x0 = [tf, u_guess, x_guess_r];
% 
% [sol, DV] = fmincon(@dv_eval, x0, [], [], [], [], lb, ub, @nonlcon, options);
% 
% 
