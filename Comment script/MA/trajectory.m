function [Y_tot,T_tot] = trajectory(launcher, site, N, orbit,options, sol)
% function to retrieve the optimal trajectory from the sol vector of
% fmincon

% same initial state of optimization
x0 = 0;
z0 = site.z0;
V0 = 1;
gamma0 = deg2rad(site.OMEGA);

% set time burning of 3 stage same as solution
launcher.st3.tb = sol(1);
% control from solution
chi = sol(2:end);

% initialize variable for cycle 
stage = fieldnames(launcher);
k = 0;
t0_stage = 0;
Y_tot = [];
T_tot = [];

% same cycle of nonlcon
for i = 1:numel(stage)
    
    m0 = launcher.(stage{i}).m0;
    Y0 = [x0 z0 V0 gamma0 m0];
    tf = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf];
    
    dt = (tf-t0_stage)/N;

for j = 1:N 
   
    t0 = (j-1)*dt + t0_stage;
    tf = j*dt + t0_stage;
    [T, Y] = ode45(@ascent_opt, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi(j+k));
    Y0 = Y(end,:);
    % assemble final state
    Y_tot = [Y_tot; Y];
    T_tot = [T_tot; T];

end
    k = j;
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    t0_stage = tf;  
end


