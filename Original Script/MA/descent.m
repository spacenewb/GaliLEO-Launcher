function [dY] = descent(~, Y, launcher, site, orbit, para)

% recalling the state
% x = Y(1);
% y = Y(2);
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);

%% ADDING WIND (supposed to be added in NED axes);

[uw, vw] = wind_generator_para(z,Y,site);

% Relative velocities (plus wind);
ur = u - uw;
vr = v - vw;
wr = w;

V_norm = norm([ur vr wr]);

%% CONSTANTS
S = para.S;                                               % [m^2]   Surface
CD = para.CD;                                             % [/] Parachute Drag Coefficient
CL = para.CL;                                            % [/] Parachute Lift Coefficient
Re = orbit.Re * 1e3;
g = orbit.g / (1 + (z/Re) )^2;                            % [N/kg] magnitude of the gravitational field at zero
m = launcher.ms;                                  % [kg] descend mass

%% ATMOSPHERE DATA

if z < 0
    z = 0;
end

[~, ~, P, rho] = atmoscoesa(z+site.z0,'none');


%% REFERENCE FRAME
% The parachutes are approximated as rectangular surfaces with the normal
% vector perpendicular to the relative velocity

t_vect = [ur vr wr];                     % Tangenzial vector
h_vect = [-vr ur 0];                     % horizontal vector    

if all(abs(h_vect) < 1e-8)
    h_vect = [-vw uw 0];
end

t_vers = t_vect/norm(t_vect);            % Tangenzial versor
h_vers = -h_vect/norm(h_vect);           % horizontal versor

n_vect = cross(t_vers, h_vers);          % Normal vector
n_vers = n_vect/norm(n_vect);            % Normal versor

if (n_vers(3) > 0)                       % If the normal vector is downward directed
    n_vect = cross(h_vers, t_vers);
    n_vers = n_vect/norm(n_vect);
end

%% FORCES

D = 0.5*rho*V_norm^2*S*CD*t_vers';       % [N] Drag vector
L = 0.5*rho*V_norm^2*S*CL*n_vers';       % [N] Lift vector
Fg = m*g*[0 0 1]';                       % [N] Gravitational Force vector
F = -D+L-Fg;                              % [N] total forces vector

%% STATE DERIVATIVES

% velocity
du = F(1)/m;
dv = F(2)/m;
dw = F(3)/m;

%% FINAL DERIVATIVE STATE ASSEMBLING

dY(1:3) = [u v w]';
dY(4) = du;
dY(5) = dv;
dY(6) = dw;

dY = dY';

%% SAVING THE QUANTITIES FOR THE PLOTS
% 
% if settings.plots
%     
%     parout.integration.t = t;
%     parout.interp.alt = -z;
%     parout.wind.body_wind = [uw, vw, ww];
%     parout.wind.NED_wind = [uw, vw, ww];
%     
%     parout.air.rho = rho;
%     parout.air.P = P;
%     
%     parout.accelerations.body_acc = [du, dv, dw];
%     
%     parout.velocities = [u, v, w];
%     
% end
