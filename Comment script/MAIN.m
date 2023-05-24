%%% MAIN %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Summary: This is the MAIN script that is responsible to call all the supporting functions and scripts to run.

% Objective: Preliminary design of a small satellite launch vehicle configuration for a 50kg payload to LEO (400 km).
%            The launcher 1st stage must be solid rocket motor and should be recoverable.

% Algorithm: The design algorithm is based on initial input values for propulsion, based on the CEA code. 
%            Optimal stating algorithm is implemented using Lagrangian Multipliers to evaluate the best 
%            possible mass distribution and staging configuration. The script then completes the sizing 
%            of the stages and propulsion systems based on empirical data of the fineness ratios. Then, mass
%            estimation is carried out using empircal formulas and data curve fitting. The script continues to 
%            iterate over the sizing and optimal staging to satisfy the fitting data within a tolerance. On 
%            reaching convergence, the configuration is finalised. The script then performs frequency analysis 
%            on the launcher design. Further, the script performs trajectory determination and optimisation based
%            on the data output of the DATCOM script. The script also evaluates the revovery predictions for the 
%            1st stageFinally, the script tries to implement control algorithm to ensure and simulate optimal
%            trajectory compliance.

% Author: GaliLEO Group

% For Launch System Design course Fall Semester 2020-2021 [Prof. Filippo Maggi]
% Department of Space Engineering, Politecnico di Milano - 18 December 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialising the script and clearing workspace
close all
clear
clc

%% LOAD DATA 
path = genpath(pwd);
addpath(path);
config 

%% STAGING, PROPELLANT BUDGET 

% initialization for the cycle
err = 1;
nmax = 50; % max number of iteration
j = 1;

% cycle for mass and staging definition 
while err > 1e-6 && j < nmax
 
% optimal staging algorithm 
[launcher,DV] = optimal_staging(data);

% Diameter and reference area 
launcher.st1.S = (data.launcher.st1.C/2)^2 * pi;
launcher.st2.C = launcher.st1.C; launcher.st3.C = launcher.st1.C;
launcher.st2.S = launcher.st1.S; launcher.st3.S = launcher.st1.S;

% propulsion data 
stage = fieldnames(launcher);
for i = 1:2
    % propulsion data from srm design
    [launcher.(stage{i}),data.prop.(stage{i})] = srm_dim(launcher.(stage{i}),data.prop.(stage{i}),data.orbit,data,i);
end
% Monopropellant dimensioning
[launcher.st3, data.prop.st3] = liquid_dim(launcher.st3, data.prop.st3, data.orbit);

% Structural fitting  
load eps_sfit % file with fit law 
eps_vect(1:2) = feval(fun_eps,[launcher.st1.mp*1e-3 launcher.st2.mp*1e-3]);  % strcutual indexes for 1 and 2 stage 
k3 = 0.6*launcher.st3.mp^(-0.15);   
eps_vect(3) = k3/(1+k3);                                                     % structural index for 3 stage

% assign eps to the struct data
data.launcher.st1.eps = eps_vect(1);
data.launcher.st2.eps = eps_vect(2);
data.launcher.st3.eps = eps_vect(3);

% evaluation of the error for the cycle
err0(1) = abs(eps_vect(1) - launcher.st1.eps);
err0(2) = abs(eps_vect(2) - launcher.st2.eps);
err0(3) = abs(eps_vect(3) - launcher.st3.eps);
err = sum(err0);

% increment cycle's counter
j = j+1;
end

%% STRUCTURE

% Lengths and Diameters
launcher.st1.L = 14.52; launcher.st2.L = 1.95; launcher.st3.L = 1.83;           % final lenght of the steps overwritten from cad
launcher.st1.L_stage = launcher.st1.L + launcher.st2.L + launcher.st3.L;
launcher.st2.L_stage = launcher.st2.L + launcher.st3.L; 
launcher.st3.L_stage = launcher.st3.L; 
launcher.st1.LoD = launcher.st1.L_stage/launcher.st1.C;                         % L/D ratio for entire rocket 

% Inertia evaluation

% Calculated as a bar element
for i = 1:3
launcher.(stage{i}).Iyyf = (1/12) * launcher.(stage{i}).m0 * launcher.(stage{i}).L_stage^2;
launcher.(stage{i}).Iyye = (1/12) * (launcher.(stage{i}).m0-launcher.(stage{i}).mp) * launcher.(stage{i}).L_stage^2;
end

% Calculated with superposition of simplified shape
[launcher] = Iyy_evaluation(launcher, data); % very similar to the rod one 

% Center of gravity evaluation
[x_cg,launcher] = CoG_evaluation(launcher, data);

% Frequency Analysis
Freq(launcher, data);


%% BASE SIMULATOR
% Ballistic uncontrolled simulator to test the dynamic
if data.simulate
    
    [traj] = simulator(data.site,launcher,data.orbit,data);
    
    traj.coeff.CD = [traj.st1.coeff.CD, traj.st2.coeff.CD, traj.st3.coeff.CD];
   
end

%% OPTIMIZATION
% evaluate the trajectory optimization

if data.optimize
    
    %- direct shooting optimization 
    [sol] = traj_opt(launcher,data);
    
    %- save data and load data to avoid optimization run every time
    % save sol.mat sol
    load sol.mat
    
    % reconstruction of optimize trajectory
    [Y_tot,T_tot] = trajectory(launcher, data.site, data.optimo.N, data.orbit, data.options, sol);
    
    % state and time vector for optimize trajectory
    Y = Y_tot;
    T = T_tot;
    
    % save data of trajectory
    %save orbit_path.mat T Y
    
end

%% COMPLETE TRAJECTORY
% trajectory of first stage reentry added and plot and data evaluation
if data.trajectory
      
      trajectory_analysis(launcher,data.orbit,data.site,data);
      
end

%% CONTROLLED SIMULLATOR 
% simulator to verify the fesibility of trajectory control with a PID on
% the TVC deflection angle
if data.control
    
    % tuner for PID 
  %  tvc_tuning(launcher)
    dt = 0.1;       % delta t interval for control 
    
    % fit gamma objective to desired delta t
    [gamma_obj] = gamma_fit(dt); 
    
    % simulator
    [traj] = control(data.site,launcher,data.orbit,data,gamma_obj,dt);
    
    % Plots
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

    figure
    plot(T,rad2deg(Y(:,4)),'b','Linewidth',2)
    title('Angles')
    xlabel('time [s]');
    ylabel('[°]');
    grid on
    hold on 
    plot(T,rad2deg(Y(:,5)),'k','Linewidth',2)
   plot(obj.T,rad2deg(obj.Y(:,4)),'r--','Linewidth',2)
    legend('\gamma','\theta','\gamma_{obj}')
    
    
end
