function [dY,parout] = ascent_opt(t, Y, launcher,site,orbit,chi)

%% define state
% x = Y(1);
z = Y(2);
V = Y(3);
gamma = Y(4);
m = Y(5);

%% assign data

Re = orbit.Re * 1e3;
g = orbit.g / (1 + (z/Re) )^2;
S = launcher.S;

%% FORCE

if t <= launcher.T_time(end) 
T = interp1(launcher.T_time, launcher.T, t);
mdot = T/(launcher.isp*orbit.g);

else 
    
    T = 0;
    mdot = 0;
    
end

% atmopshere data
if z<0
    z=0;
end

[~, a, P, rho] = atmoscoesa(site.z0 + z,'none');
Mach = V/a;

CD = 0.49;
alpha = 0;   
if isnan(rho)
    rho = 0;
end
    
% aerodynamics force
q_dyn = 0.5*rho*V^2;            %dynamic pressure
D = q_dyn*S*CD;
L = 0;


%% ASSEMBLE DERIVATES
r = z+Re;

dV = (T*cos(chi) -D)/m -g*sin(gamma);
% if V < 10 
%     dgamma = 0;
% else
    dgamma = (V*cos(gamma))/r + (T*sin(chi) + L)/(m*V) - (g*cos(gamma))/V;
% end

dz = V*sin(gamma);

dY(1) = V*cos(gamma);
dY(2) = dz;
dY(3) = dV;
dY(4) = dgamma;
dY(5) = - mdot;

dY = dY';

%% Parout 

parout.interp.M = Mach;
parout.interp.alpha = alpha;

parout.forces.AeroDyn_Forces = [D, L, T];
parout.accelerations.body_acc = dV;

parout.air.rho = rho;
parout.air.P = P;

parout.coeff.CD = CD;
