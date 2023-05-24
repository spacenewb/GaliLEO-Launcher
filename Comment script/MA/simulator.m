function [traj] = simulator(site,launcher,orbit,data)

options = data.options;
options_ascent = data.options_ascent;
options_descent = data.options_descent;



%% Initial condition

x0 = 0;
z0 = site.z0;
V0 = 34.9;
gamma0 = deg2rad(site.OMEGA);
omega0 = - site.lon0;


%% Ode solver

% initialize variable for cycle 
stage = fieldnames(launcher);
Y_tot = [];
T_tot = [];
t0 = 0;
chi = 0;

for i = 1:numel(stage)
    
    m0 = launcher.(stage{i}).m0;
    Y0 = [x0 z0 V0 gamma0 m0 omega0];
    tf = launcher.(stage{i}).tb + t0;
    launcher.(stage{i}).T_time = [t0 tf];
    
  
    [T, Y] = ode45(@rocket_dyn, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi);
    
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    omega0 = Y(end,6);
    
    Y_tot = [Y_tot; Y];
    T_tot = [T_tot; T];
    t0 = tf;
    
    [traj.(stage{i})] = RecallOde(@rocket_dyn, T, Y, launcher.(stage{i}), site, orbit, chi);
    
end

%% Assembly final state 

traj.T = T_tot;
traj.Y = Y_tot;



%% First stage dynamics 

% % ascent 
% launcher.st1.T_time = [0 launcher.st1.tb];
% Y01 = Y_tot(T_tot==launcher.st1.tb,:)';
% Y01 = Y01(:,1);
% Y01(5) = launcher.st1.ms;
% [T1a, Y1a] = ode45(@rocket_dyn, [launcher.st1.tb 2000], Y01, options_ascent, launcher.st1, site, orbit, chi);
% [traj.bal] = RecallOde(@rocket_dyn, T1a, Y1a, launcher.st1, site, orbit, chi);
% 
% % descent 
% t0 = T1a(end);
% x0d = Y1a(end,1)*cosd(site.azimuth);
% y0d = Y1a(end,1)*sind(site.azimuth);
% z0d = Y1a(end,2);
% vxd = Y1a(end,3)*cosd(site.azimuth)*cos(Y1a(end,4));
% vyd = Y1a(end,3)*sind(site.azimuth)*cos(Y1a(end,4));
% vzd = Y1a(end,3)*sin(Y1a(end,4));
% T1d_tot = [];
% Y1d_tot = [];
% 
% for i = 1 : length(data.para)
%     
% Y01d = [x0d y0d z0d vxd vyd vzd];
% [T1d, Y1d] = ode45(@descent, [t0 6000], Y01d, options_descent, launcher.st1, site, orbit, data.para(i));
% 
% x0d = Y1d(end,1);
% y0d = Y1d(end,2);
% z0d = Y1d(end,3);
% vxd= Y1d(end,4);
% vyd = Y1d(end,5);
% vzd = Y1d(end,6);
% t0 = T1d(end);
% 
% T1d_tot = [T1d_tot; T1d];
% Y1d_tot = [Y1d_tot; Y1d];
% 
% end


% traj.T1a = 0; % T1a;
% traj.Y1a =  0; % Y1a;
% traj.T1d =  0; %T1d_tot;
% traj.Y1d =  0; %Y1d_tot;


