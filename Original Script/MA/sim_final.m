function [traj] = sim_final(site,launcher,orbit,options,gamma_obj,dt)

x0 = 0;
z0 = site.z0;
V0 = 5;
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
Kp = - 2.1;
Kd = - 0.7;
Ki =  - 0.09;
%%%%%%%%%%%%%%%%%%%%
% Kp1 = - 6;
% Kp2 = - 2;
% Kd = - 0.5;
% Ki1 =  - 0.8;
% Ki2 =  Ki1;



tic
% i = 1;
Y_tot = [];
T_tot = [];

theta = gamma_obj(1);
gamma_old = gamma_obj(1);
int_e = 0;
e_old = 0;

for i = 1:3
    
    m0 = launcher.(stage{i}).m0;
    Iyy0 = launcher.(stage{i}).Iyyf;
    Y0 = [x0 z0 V0 gamma0 theta0 w0 m0 Iyy0 omega0];
    tf_stage = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf_stage];
    
    if tf_stage > (length(gamma_obj)*dt)
        tf_stage = (length(gamma_obj)*dt);
    end
    

% tf_stage = 10; 
% Y_tot = zeros(length(t0_stage:dt:tf_stage-dt),9);
% T_tot = zeros(length(t0_stage:dt:tf_stage-dt),1);
 for j = t0_stage:dt:tf_stage-dt
    
     t0 = j;
     tf = j+dt;
          %%%% PID
     
    e =  gamma_obj(k) - theta;
    dot_gamma = (gamma_obj(k) - gamma_old)/dt;
   dot_e = dot_gamma - Y0(6);
    int_e = int_e+(e+e_old)*dt/2;
    chi = Kp*e + Kd*dot_e + Ki*int_e;
    
 %  chi = Kp * (e+Td*dot_e + int_e/Ti)
     if abs(chi) > deg2rad(7)
         chi = deg2rad(7)*(chi/abs(chi));
     end

    [T, Y] = ode45(@rocket_fin, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi);
    
    theta = Y(end,5);
    Y0 = Y(end,:);
    gamma_old = gamma_obj(k);
%    e_old_gamma = e_gamma;
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

