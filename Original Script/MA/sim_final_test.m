function [traj] = sim_final_test(site,launcher,orbit,options,gamma_obj,dt)

x0 = 0;
z0 = site.z0;
vx0 = 0;
vz0 = 0;
theta0 =  deg2rad(site.OMEGA);
p0 = 0;

sol = round(length(gamma_obj)*dt - (launcher.st1.tb + launcher.st2.tb));
launcher.st3.tb = sol;

stage = fieldnames(launcher);
t0_stage = 0;
k = 1;

% stage 1 e 2
Kp = - 3;
Kd = -0.25;
Ki =  - 0.09;
% stage 3


T_tot = [];
Y_tot = [];
theta = gamma_obj(1);
gamma_old = gamma_obj(1);
int_e = 0;
e_old = 0;


tic
i = 3;
 for i = 1:2
    
    m0 = launcher.(stage{i}).m0;
    Iyy0 = launcher.(stage{i}).Iyyf;
    Y0 = [x0 z0 vx0 vz0 theta0 p0 m0 Iyy0];
    tf_stage = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf_stage];
% if tf_stage > (length(gamma_obj)*dt)
%     tf_stage = (length(gamma_obj)*dt);
% end
% load sate2.mat
% Kp = -0.00914; Ki = -0.000529; Kd = -0.0395;
% tf_stage = t0_stage + 5;
% theta = Y0(5);
% gamma_old = Y0(5);
% int_e = 0;
% e_old = 0;

% if i == 3
%     save sate2.mat
% end

% Y_tot = zeros(length(t0_stage:dt:tf_stage-dt),8);
% T_tot = zeros(length(t0_stage:dt:tf_stage-dt),1);
 for j = t0_stage:dt:tf_stage-dt
    
     t0 = j;
     tf = j+dt;
    
    %%% PID
    
    e = gamma_obj(k) - theta;
    dot_gamma = gamma_obj(k) - gamma_old;
    dot_e = dot_gamma - Y0(6);
    int_e = int_e+(e+e_old)*dt/2;
    chi = Kp*e + Ki*int_e + Kd*dot_e;
    
         if abs(chi) > deg2rad(7)
             chi = deg2rad(7)*(chi/abs(chi));
         end

    [T, Y] = ode45(@rocket_fin_test, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi);
    
    theta = Y(end,5);
    gamma_old = gamma_obj(k);
    Y0 = Y(end,:);
%     Y_tot(k,:) = Y0;
%     T_tot(k,1) = T(end);
    e_old = e;
    k = k +1;
    T_tot = [T_tot; T];
    Y_tot = [Y_tot;Y];
    
 end


    x0 = Y(end,1);
    z0 = Y(end,2);
    vx0 = Y(end,3);
    vz0 = Y(end,4);
    theta0 = Y(end,5);
    p0 = Y(end,6);
    t0_stage = tf;
 
   
  end
toc
traj.Y = Y_tot;
traj.T = T_tot;

