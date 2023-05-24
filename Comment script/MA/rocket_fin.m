function [dY, dgamma, alpha] = rocket_fin(t, Y, launcher,site,orbit,delta)

% x = Y(1);
z = Y(2);
V = Y(3);
gamma = Y(4);
theta = Y(5);
w = Y(6);
m = Y(7);
Iyy = Y(8);

%% assign data

Re = orbit.Re * 1e3;
g = orbit.g / (1 + (z/Re) )^2;
S = launcher.S;
C = launcher.C;
tb = launcher.tb;
L_stage = launcher.L_stage;


% Adding wind
[uw, vw] = wind_generator(z,Y,site);
wind = norm([uw vw]);

% relative velocities in intertial frame
u = V*cos(gamma);
v = V*sin(gamma);

ur = u - wind;



% evaluate angle of attack
if not(u < 1e-9 ||  V < 1e-9)
    alpha = theta - atan2(v,ur);
else
    alpha = 0;
end


%% FORCE

if t <= launcher.T_time(end) 
T = interp1(launcher.T_time, launcher.T, t);
mdot = T/(launcher.isp*orbit.g);
Idot = (launcher.Iyyf - launcher.Iyye)/tb;
xcg = interp1(launcher.T_time,[launcher.xcgf launcher.xcge], t);    % xcg taken from the tip of launcher


else 
    
    T = 0;
    mdot = 0;
    Idot = 0;
    xcg = launcher.xcge;
    
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
CMAfull = launcher.aero.full.Coeffs.CMA;


CDempty = launcher.aero.empty.Coeffs.CD;
CLempty = launcher.aero.empty.Coeffs.CL;
CMAempty = launcher.aero.empty.Coeffs.CMA;
    
% For interpolation angles are in degree!
    
givA = launcher.aero.full.State.Alphas*pi/180;
givB = launcher.aero.full.State.Betas*pi/180;
givAlts = launcher.aero.full.State.Altitudes;
givM = launcher.aero.full.State.Machs;

% INTERPOLATION AT THE BOUNDARIES

z_d = z;
alpha_d = alpha;
if Mach > givM(end)
    
    Mach = givM(end);
    
end

if Mach < givM(1)
    
    Mach = givM(1);
    
end

if alpha > givA(end)
    
    alpha_d = givA(end);
    
elseif alpha < givA(1)
    
    alpha_d = givA(1);
    
end


if z > givAlts(end)
    
    z_d = givAlts(end);
    
elseif z < givAlts(1)
    
    z_d = givAlts(1);
    
end

    
    %Multi-variable interpolation (with 'spline' option interpn extrapolates
    %also outside the boundaries of the data set);
    
    CDf=interp4_easy(givA,givM,givB,givAlts,CDfull,alpha_d,Mach,0,z_d);
    CLf=interp4_easy(givA,givM,givB,givAlts,CLfull,alpha_d,Mach,0,z_d);
    CMAf=interp4_easy(givA,givM,givB,givAlts,CMAfull,alpha_d,Mach,0,z_d);
    
    CDe=interp4_easy(givA,givM,givB,givAlts,CDempty,alpha_d,Mach,0,z_d);
    CLe=interp4_easy(givA,givM,givB,givAlts,CLempty,alpha_d,Mach,0,z_d);
    CMAe=interp4_easy(givA,givM,givB,givAlts,CMAempty,alpha_d,Mach,0,z_d);

    
    %Linear variation from full to empty configuration
    if t < launcher.T_time(end) 
        CD = (t-launcher.T_time(1))/tb*(CDe-CDf)+CDf;
        CL = (t-launcher.T_time(1))/tb*(CLe-CLf)+CLf;
        CMA = (t-launcher.T_time(1))/tb*(CMAe-CMAf)+CMAf;
        

    else
        CD = CDe;
        CL = CLe;
        CMA = CMAe;
        
    end
    % aerodynamics force
    q_dyn = 0.5*rho*V^2;            %dynamic pressure
    D = q_dyn*S*CD;
    L = q_dyn*S*CL;
    M = q_dyn*S*C*CMA*alpha;
    
else
    
    CD = NaN;
    D = 0;
    L = 0;
    M = 0;
    
end

%% ASSEMBLE DERIVATES
r = z+Re;
% gamma = theta+alpha;

dV = (T*cos(alpha-delta) -D)/m -g*sin(gamma);
dgamma = (V*cos(gamma))/r + (T*sin(alpha-delta) + L)/(m*V) - (g*cos(gamma))/V;
dw = (M + T*sin(delta)*(xcg-L_stage))/Iyy; 


dz = V*sin(gamma);
domega = (V*cos(gamma))/(r);

dY(1) = V*cos(gamma);
dY(2) = dz;
dY(3) = dV;
dY(4) = dgamma;
dY(5) = w;
dY(6) = dw;
dY(7) = - mdot;
dY(8) = - Idot;
dY(9) = - domega;


dY = dY';
