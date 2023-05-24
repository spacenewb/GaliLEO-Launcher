%%% MAX CROSSWIND EVALUATOR %%%


run MAIN 
%tvc_tuning(launcher)

Y0(1) = 0;
Y0(2) = 0;
Y0(3) = 0;
Y0(4) = 0;
Y0(5) = deg2rad(90);
Y0(6) = 0;

tspan = [0 20];
wind = -30; %[m/s]
s = length(wind);
chi = deg2rad(0);
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);
Tt = [];
Yt = [];

Kp = -15;
Ki = -0.9623;
Kd = -2.66;

dt = 0.1;
theta = deg2rad(90);
p = 0;
dot_gamma = 0;
int_e = 0;
e_old = 0;
k = 1;
for t = 0:dt:(2-dt)
    
    
    
    % phis-phi;
    e = deg2rad(90) - theta;
    
    dot_gamma = 0;
    dot_e = dot_gamma - p;
    int_e = int_e+(e+e_old)*dt/2;
    
%     delta = kp*(e+Td*dot_e+int_e/Ti);
      chi = Kp*e + Kd*dot_e + Ki*int_e;
      
     if abs(chi) > deg2rad(7)
         chi = deg2rad(7)*(chi/abs(chi));
     end
     
     chiVect(k) = rad2deg(chi);
    
    [T,Y] = ode45(@start_dynamic, [t t+dt], Y0, options, launcher, chi, wind);
    Y0 = Y(end,:);
    
    
    theta = Y(end,5);
    p = Y(end,6);
    e_old = e;
    Yt = [Yt;Y];
    Tt = [Tt;T];
    
k = k +1;
end


figure()
plot(Tt,rad2deg(Yt(:,5)))

figure()
plot(Tt,Yt(:,1))
title('time vs x')

figure()
plot(Yt(:,1),Yt(:,2))
hold on
plot([-5 -5],[0 20],'k--','LineWidth',5) 
xlabel('x')
ylabel('z')















% delta_max = deg2rad(6.5); 
% sin_delta = sin(delta_max); % sin(delta) [rad]
% 
% L_tot = launcher.st1.L_stage; %Lb [m]
% 
% x_cg = launcher.st1.xcgf; % x_cg [m from the tip]
% 
% x_cp = x_cg - 5.5*0.83; % x_cp [m from the tip], x_cp set as less than xcg to be into unstable conditions (worst ones)
% 
% S = launcher.st1.S; %S [m^2]
% 
% T = launcher.st1.T(1); % T [N]
% 
% Cn_alpha = 3.3; % --> DATA COMING FROM AERODYNAMICS
% 
% v = [24:1:50]'; %vector of velocities [m/s]
% rho = 1.225; % [kg/m^3] assumed as the density of the air at z = 0m
% q = 0.5*rho*(v.^2); % q 
% 
% theta_20m = rad2deg(Y(end,5));
% sin_theta = sin(Y(end,5));
% 
% k = length(v);
% 
% for i = 1:k
%     
%     alpha_w(i) = -(sin_delta*T*(L_tot-x_cg))/(Cn_alpha*S*q(i)*(x_cg-x_cp));%[rad]
%     sin_alpha_w(i) = sin(alpha_w(i));
%     %v_wind (i) = v(i)*(sin_alpha_w(i)/sin_theta);
%     v_wind (i) = atan(alpha_w(i))*v(i);
% end
% 
% max_alpha_w = rad2deg(max(abs(alpha_w)));
% max_v_wind = min(abs(v_wind));
% 
% 
% 
% 
% 
