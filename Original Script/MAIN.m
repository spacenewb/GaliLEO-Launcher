%%% MAIN %%%
close all
clear
clc

%% LOAD DATA 
path = genpath(pwd);
addpath(path);
config 

%% STAGING, PROPELLANT BUDGET AND STRUCTURE
err = 1;
nmax = 50;
j = 1;

while err > 1e-6 && j < nmax
    

[launcher,DV] = optimal_staging(data);

% Diameter and reference area 
launcher.st1.S = (data.launcher.st1.C/2)^2 * pi;
launcher.st2.C = launcher.st1.C; launcher.st3.C = launcher.st1.C;
launcher.st2.S = launcher.st1.S; launcher.st3.S = launcher.st1.S;

% propulsion data 
stage = fieldnames(launcher);
for i = 1:2
    
    [launcher.(stage{i}),data.prop.(stage{i})] = srm_dim(launcher.(stage{i}),data.prop.(stage{i}),data.orbit,data,i);
     
end
[launcher.st3, data.prop.st3] = liquid_dim(launcher.st3, data.prop.st3, data.orbit);

% Structural part 
load eps_sfit
% [ms_vect, eps_vect] = structural_masses (launcher, data.prop, data);
eps_vect(1:2) = feval(fun_eps,[launcher.st1.mp*1e-3 launcher.st2.mp*1e-3]);
% eps_vect(1:2) = [5.3366 * launcher.st1.mp^(-0.4146), 5.3366 * launcher.st2.mp^(-0.4146)];
k3 = 0.6*launcher.st3.mp^(-0.15);
eps_vect(3) = k3/(1+k3);

data.launcher.st1.eps = eps_vect(1);
data.launcher.st2.eps = eps_vect(2);
data.launcher.st3.eps = eps_vect(3);

%err = abs(max(eps_vect - [launcher.st1.eps launcher.st2.eps launcher.st3.eps]));
err0(1) = abs(eps_vect(1) - launcher.st1.eps);
err0(2) = abs(eps_vect(2) - launcher.st2.eps);
err0(3) = abs(eps_vect(3) - launcher.st3.eps);
err = sum(err0);
j = j+1;

% hold on
% plot(j,err,'o-')

end

% Lengths and Diameters
launcher.st1.S = (data.launcher.st1.C/2)^2 * pi;
launcher.st2.C = launcher.st1.C; launcher.st3.C = launcher.st1.C;
launcher.st2.S = launcher.st1.S; launcher.st3.S = launcher.st1.S;
launcher.st1.L = 14.52; launcher.st2.L = 1.95; launcher.st3.L = 1.83; % from cad 
launcher.st1.L_stage = launcher.st1.L + launcher.st2.L + launcher.st3.L;
launcher.st2.L_stage = launcher.st2.L + launcher.st3.L; 
launcher.st3.L_stage = launcher.st3.L; 
launcher.st1.LoD = launcher.st1.L_stage/launcher.st1.C;

% Inertia 

% preliminary could be calculated as: 
for i = 1:3
launcher.(stage{i}).Iyyf = (1/12) * launcher.(stage{i}).m0 * launcher.(stage{i}).L_stage^2;
launcher.(stage{i}).Iyye = (1/12) * (launcher.(stage{i}).m0-launcher.(stage{i}).mp) * launcher.(stage{i}).L_stage^2;
end

%[launcher] = Iyy_evaluation(launcher, data); % very similar to the rod one 

% Center of gravity evaluation

%[x_cg,launcher] = CoG_evaluation(launcher, data);

% Frequency Analysis5
%Freq(launcher, data);


%% BASE SIMULATOR

if data.simulate
    
    [traj] = simulator(data.site,launcher,data.orbit,data);
    
    traj.air.P = [traj.st1.air.P traj.st2.air.P traj.st3.air.P];
    traj.air.rho = [traj.st1.air.rho traj.st2.air.rho traj.st3.air.rho];
    traj.accelerations = [traj.st1.accelerations.body_acc, traj.st2.accelerations.body_acc, traj.st3.accelerations.body_acc];
    traj.interp.M = [traj.st1.interp.M traj.st2.interp.M traj.st3.interp.M];
    traj.interp.alpha = [traj.st1.interp.alpha traj.st2.interp.alpha traj.st3.interp.alpha];
    traj.aero_force = [traj.st1.forces.AeroDyn_Forces; traj.st2.forces.AeroDyn_Forces; traj.st3.forces.AeroDyn_Forces];
    traj.coeff.CD = [traj.st1.coeff.CD, traj.st2.coeff.CD, traj.st3.coeff.CD];
    
    % PLOT TRAJECTORY
    if data.plot_traj
        plot_sim(traj,data);
    end
    
    % evaluate loss 
    Y = traj.Y;
    T = traj.T;
    Re = data.orbit.Re * 1e3;
    drag_loss = traj.aero_force(:,1)./Y(:,5);
    dv_drag = trapz(T,drag_loss);
    grav_loss = data.orbit.g .* ( Re./(Re + Y(:,2)) ).^2 .* sin(Y(:,4));
    dv_grav = trapz(T,grav_loss);
    
    dv_loss = dv_grav+dv_drag;
    
    
    
end

%% OPTIMIZATION 

if data.optimize
    
   % [sol] = traj_opt(launcher,data);
    N = data.optimo.N;
  %  save sol.mat sol
  load sol.mat 
    [Y_tot,T_tot] = trajectory(launcher, data.site, N, data.orbit, data.options, sol);
    
    Y = Y_tot;
    T = T_tot;
    
    figure
    plot(T,Y(:,2)/1000,'b','Linewidth',2)
    title('Altitude Profile vs Time')
    xlabel('time [s]');
    ylabel('z [km]');
    grid on
    
    figure
    plot(Y(:,1)/1000,Y(:,2)/1000,'b','Linewidth',2)
    title('Trajectory')
    xlabel('x [km]');
    ylabel('z [km]');
    grid on 
    
    
     figure
    plot(T,Y(:,3).*cos(Y(:,4)),'b','Linewidth',2)
    title('vx')
    xlabel('time [s]');
    ylabel('vx [m/s]');
    grid on
    
    figure
    plot(T,Y(:,3).*sin(Y(:,4)),'b','Linewidth',2)
    title('vz')
    xlabel('time [s]');
    ylabel('vz [m/s]');
    grid on
    
    figure
    plot(T,Y(:,3),'b','Linewidth',2)
    title('V')
    xlabel('time [s]');
    ylabel('V [m/s]');
    grid on
   
       
    figure
    plot(T,rad2deg(Y(:,4)),'b','Linewidth',2)
    title('Attitude')
    xlabel('time [s]');
    ylabel('theta [°]');
    grid on
    
%     save orbit_path.mat T Y
    
end

%% FIRST STAGE SIMULATOR

if data.recovery
      
      recovery_sim(launcher,data.orbit,data.site,data)
      
end



%% CONTROLLED SIMULLATOR 

if data.simulate_final
    
  %  tvc_tuning(launcher)
    dt = 0.1;
    [gamma_obj] = gamma_fit(dt);
    [traj] = sim_final(data.site,launcher,data.orbit,data,gamma_obj,dt);
    
obj = load('orbit_path.mat');
    Y = traj.Y;
    T = traj.T;
    
    figure
    plot(Y(:,1)/1000,Y(:,2)/1000,'k','Linewidth',2)
    title('Altitude vs Range')
    xlabel('Range [km]');
    ylabel('Altitude [km]');
    grid on
    hold on 
    plot(obj.Y(:,1)/1000,obj.Y(:,2)/1000,'r--','Linewidth',2)
    legend('Controlled Trajectory','Objective trajectory')
%     
%      figure
%     plot(T,sqrt(Y(:,3).^2 + Y(:,4).^2),'b','Linewidth',2)
%     title('vx')
%     xlabel('time [s]');
%     ylabel('vx [m/s]');
%     grid on
%     
%     
%     figure
%     plot(T,Y(:,3),'b','Linewidth',2)
%     title('vz')
%     xlabel('time [s]');
%     ylabel('vz [m/s]');
%     grid on
%     
%     figure
%     plot(T,Y(:,4),'b','Linewidth',2)
%     title('vz')
%     xlabel('time [s]');
%     ylabel('vz [m/s]');
%     grid on
%     
%     figure
%     plot(T,Y(:,6),'b','Linewidth',2)
%     title('p')
%     xlabel('time [s]');
%     ylabel('p [m/s]');
%     grid on
%     
    figure
    plot(T,rad2deg(Y(:,4)),'b','Linewidth',2)
    title('Angles')
    xlabel('time [s]');
    ylabel('[°]');
    grid on
    hold on 
    plot(T,rad2deg(Y(:,5)),'k','Linewidth',2)
   % plot(0:dt:obj.T(end),rad2deg(gamma_obj),'g--','Linewidth',2)
   plot(obj.T,rad2deg(obj.Y(:,4)),'r--','Linewidth',2)
    legend('\gamma','\theta','\gamma_{obj}')
    
    
end


%clearvars data
