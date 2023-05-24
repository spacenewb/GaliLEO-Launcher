%%% MAX CROSSWIND EVALUATOR %%%
%
% ASSUMPTION: gamma = theta, so alpha = alpha_wind
%             For the take-off transitory, set as a 20m meter ascent phase,
%             the control angle is assumed to be null and the rocket "stabilized" only by its big inertia. 
%             Then afer the first
%             20m the control starts to apply from the angle theta obtained
%             after the first 20m.
%
%if we have a wind disturbance alpha_w, the TVC (angle delta) has a minimun deflection
%required:
%
% sin(delta) = - (Cn_alha*alpha_w*S*q*(x_cp - x_cg))/(T*(Lb - x_cg)) ;
%
% I have to finde the maximum alpha_w for which delta remains into its
% boundary values [-6.5° ; +6.5°]
%
% DATA NEEDED:
%
% Cn_alpha = Normal force coefficent by the angle of attack 
% alpha_w = under hour hP is the same as the angle of attack --> OUTPUT
% S = cross surface
% q = dynamic pressure
% x_cg = baricenter position
% x_cp = center of pressure position
% Lb = length of the body (whole)
% T = thrust
% sin(delta) = sin of the maximum deflection allowable from TVC
%
run MAIN
Y0(1) = 0;
Y0(2) = 0;
Y0(3) = 0;
Y0(4) = 0;
Y0(5) = deg2rad(90);
Y0(6) = 0;

tspan = [0:1:30];
wind = -40; %[m/s] first guess on wind speed 
s = length(wind);
chi = deg2rad(0);%setting of the control angle as zero for the first ascent transitory
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8,'Events',@stop);


winds = wind;
[T,Y] = ode45(@start_dynamic, tspan, Y0, options, launcher, chi, winds);%dynamic for the first ascent transitory

delta_max = deg2rad(6.5); 
sin_delta = sin(delta_max); % sin(delta) [rad]

L_tot = launcher.st1.L_stage; %Lb [m]

x_cg = launcher.st1.xcgf; % x_cg [m from the tip]

x_cp = x_cg - 5.5*0.83; % x_cp [m from the tip], x_cp set as less than xcg to be into unstable conditions (worst ones)

S = launcher.st1.S; %S [m^2]

T = launcher.st1.T(1); % T [N]

Cn_alpha = 3.3; % --> DATA COMING FROM AERODYNAMICS

v = [24:1:50]'; %vector of velocities [m/s]
rho = 1.225; % [kg/m^3] assumed as the density of the air at z = 0m
q = 0.5*rho*(v.^2); % q 

theta_20m = rad2deg(Y(end,5)); %theta angle after the 20m ascent transitory without control
sin_theta = sin(Y(end,5));

k = length(v);

%maximum wind speed admissibile after 20m ascent:
for i = 1:k
    
    alpha_w(i) = -(sin_delta*T*(L_tot-x_cg))/(Cn_alpha*S*q(i)*(x_cg-x_cp));%[rad]
    sin_alpha_w(i) = sin(alpha_w(i));
    v_wind (i) = atan(alpha_w(i))*v(i);
end

max_alpha_w = rad2deg(max(abs(alpha_w))); %maximum wind angle admissible
max_v_wind = min(abs(v_wind)); %maximum admissible wind speed at take off