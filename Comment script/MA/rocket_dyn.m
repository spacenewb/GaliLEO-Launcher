function [dY, parout] = rocket_dyn(t, Y, launcher,site,orbit,chi)

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
tb = launcher.tb ;

% Adding wind
[uw, vw] = wind_generator(z,Y,site);
wind = norm([uw vw]);

% relative velocities in body axes
u = V*cos(gamma);
v = V*sin(gamma);

ur = u - wind(1);

% evaluate angle of attack
if not(u < 1e-9 ||  V < 1e-9)
    alpha = gamma -  atan2(v,ur);
else
    alpha = 0;
end


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


if not(isnan(rho)) 
    
% assign aerodynamic coefficient

CDfull = launcher.aero.full.Coeffs.CD;
CLfull = launcher.aero.full.Coeffs.CL;

CDempty = launcher.aero.empty.Coeffs.CD;
CLempty = launcher.aero.empty.Coeffs.CL;
    
% For interpolation angles are in degree!
    
givA = launcher.aero.full.State.Alphas*pi/180;
givB = launcher.aero.full.State.Betas*pi/180;
givAlts = launcher.aero.full.State.Altitudes;
givM = launcher.aero.full.State.Machs;

% INTERPOLATION AT THE BOUNDARIES

z_d = z;
if Mach > givM(end)
    
    Mach = givM(end);
    
end

if Mach < givM(1)
    
    Mach = givM(1);
    
end

if alpha > givA(end)
    
    alpha = givA(end);
    
elseif alpha < givA(1)
    
    alpha = givA(1);
    
end


if z > givAlts(end)
    
    z_d = givAlts(end);
    
elseif z < givAlts(1)
    
    z_d = givAlts(1);
    
end

    
    %Multi-variable interpolation (with 'spline' option interpn extrapolates
    %also outside the boundaries of the data set);
    
    CDf=interp4_easy(givA,givM,givB,givAlts,CDfull,alpha,Mach,0,z_d);
    CLf=interp4_easy(givA,givM,givB,givAlts,CLfull,alpha,Mach,0,z_d);
    
    CDe=interp4_easy(givA,givM,givB,givAlts,CDempty,alpha,Mach,0,z_d);
    CLe=interp4_easy(givA,givM,givB,givAlts,CLempty,alpha,Mach,0,z_d);

    
    %Linear variation from full to empty configuration
    if t < launcher.T_time(end) 
        CD = t/tb*(CDe-CDf)+CDf;
        CL = t/tb*(CLe-CLf)+CLf;

    else
        CD = CDe;
        CL = CLe;
        
    end
    % aerodynamics force
    q_dyn = 0.5*rho*V^2;            %dynamic pressure
    D = q_dyn*S*CD;
    L = q_dyn*S*CL;
    
else
    
    CD = NaN;
    D = 0;
    L = 0;
    
end

%% ASSEMBLE DERIVATES
r = z+Re;

dV = (T*cos(chi) -D)/m -g*sin(gamma);
dgamma = (V*cos(gamma))/r + (T*sin(chi) + L)/(m*V) - (g*cos(gamma))/V;
dz = V*sin(gamma);
domega = (V*cos(gamma))/(r);

dY(1) = V*cos(gamma);
dY(2) = dz;
dY(3) = dV;
dY(4) = dgamma;
dY(5) = - mdot;
dY(6) = - domega;

dY = dY';

%% SAVING THE QUANTITIES FOR THE PLOTS


parout.interp.M = Mach;
parout.interp.alpha = alpha;

parout.forces.AeroDyn_Forces = [D, L, T];
parout.accelerations.body_acc = dV;

parout.air.rho = rho;
parout.air.P = P;

parout.coeff.CD = CD;
% parout.coeff.XCP = XCP_value;





