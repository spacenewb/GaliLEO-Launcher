function [traj] = sim_final_gamma(site,launcher,orbit,options,gamma_obj,dt)

x0 = 0;
z0 = site.z0;
V0 = 1;
gamma0 = deg2rad(site.OMEGA);
theta0 = gamma0;
w0 = 0;
omega0 = - site.lon0;

sol = round(length(gamma_obj)*dt - (launcher.st1.tb + launcher.st2.tb));
launcher.st3.tb = sol;

stage = fieldnames(launcher);
t0_stage = 0;
k = 1;

%%%%% TOP %%%%%
% Kp = - 2.1;
% Kd = - 0.7;
% Ki =  - 0.09;
%%%%%%%%%%%%%%%%%%%%
Kp1 = - 6;
Kp2 = - 7;
Kd1 = - 2;
Kd2 = - 6;
Ki1 =  - 0.9 ;
Ki2 =  - 11.5;

% Kp1 =  -1;
% Kp2 =  1;
% Kd1 = - 1.5;
% Kd2 =  1.5;
% Ki1 =  - 0.5 ;
% Ki2 =   0.5;




tic
% i = 1;
Y_tot = [];
T_tot = [];

theta = gamma_obj(1);
gamma_old = gamma_obj(1);
gamma = gamma_obj(1);
gamma_prev = gamma_obj(1);
int_e_t = 0;
int_e_g = 0;
e_old_theta = 0;
e_old_gamma = 0;

for i = 1:1
    
    m0 = launcher.(stage{i}).m0;
    Iyy0 = launcher.(stage{i}).Iyyf;
    Y0 = [x0 z0 V0 gamma0 theta0 w0 m0 Iyy0 omega0];
    tf_stage = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf_stage];
    
    if tf_stage > (length(gamma_obj)*dt)
        tf_stage = (length(gamma_obj)*dt);
    end
    

      tf_stage = 10; 
% Y_tot = zeros(length(t0_stage:dt:tf_stage-dt),9);
% T_tot = zeros(length(t0_stage:dt:tf_stage-dt),1);
 for j = t0_stage:dt:tf_stage-dt
    
     t0 = j;
     tf = j+dt;
         %%%% PID
    
    e_theta =  gamma - theta;
    e_gamma =  gamma_obj(k) - gamma;
    dot_gamma_obj = (gamma_obj(k) - gamma_old)/dt;
    dot_gamma = (gamma - gamma_prev)/dt;
    dot_e_theta = dot_gamma_obj - Y0(6);
    dot_e_gamma = dot_gamma_obj - dot_gamma;
    int_e_t = int_e_t+(e_theta+e_old_theta)*dt/2;
    int_e_g = int_e_g+(e_gamma+e_old_gamma)*dt/2;
    delta1 = Kp1*e_theta + Kd1*dot_e_theta + Ki1*int_e_t;
    delta2 = Kp2*e_gamma + Kd2*dot_e_gamma + Ki2*int_e_g;
    chi = (delta1+delta2)/2;
    
     if abs(chi) > deg2rad(7)
         chi = deg2rad(7)*(chi/abs(chi));
     end

    [T, Y] = ode45(@rocket_fin, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi);
    
    
    gamma_prev = gamma;
    theta = Y(end,5);
    gamma = Y(end,4);
    gamma_old = gamma_obj(k);
    Y0 = Y(end,:);
    e_old_theta = e_theta;
    e_old_gamma = e_gamma;
    k = k + 1;

T_tot = [T_tot; T];
Y_tot = [Y_tot;Y];
 end


    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    theta0 = Y(end,5);
    w0 = Y(end,6);
    omega0 = Y(end,9);
    t0_stage = tf;
 
   
end
toc
traj.Y = Y_tot;
traj.T = T_tot;
