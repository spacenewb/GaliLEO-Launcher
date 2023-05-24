function [dY, parout] = rocket_fin_test(t, Y, launcher,site,orbit,chi)

% x = Y(1);
z = Y(2);
vx = Y(3);
vz = Y(4);
theta = Y(5);
p = Y(6);
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

vxr = vx - wind;
V = norm(vxr,vz);

% evaluate angle of attack
if not(vx < 1e-9 ||  V < 1e-9)  
    gamma = atan(vz/vxr);
    alpha = theta - gamma;
else
    gamma = theta;
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
    A = q_dyn*S*CD;
    N = q_dyn*S*CL;
    M = q_dyn*S*C*CMA*alpha;
    
else
    
    CD = NaN;
    A = 0;
    N = 0;
    M = 0;
    
end

%% ASSEMBLE DERIVATES
r = z+Re;
Fx = T*cos(theta+chi)-A*cos(theta)-N*sin(theta); %- (V^2*cos(gamma)/r) * sin(gamma);
Fy = T*sin(theta+chi)-A*sin(theta)+N*cos(theta)-m*g + (V^2*cos(theta)/r);
dp = (M + T*sin(chi)*(xcg-L_stage))/Iyy;

dY(1) = vx;
dY(2) = vz;
dY(3) = (Fx)/m ;
dY(4) = (Fy)/m;
dY(5) = p;
dY(6) = dp;
dY(7) = - mdot;
dY(8) = - Idot;

dY = dY';

if isnan(dY)
    pause
end

%% SAVING THE QUANTITIES FOR THE PLOTS

% 
% parout.interp.M = Mach;
% parout.interp.alpha = alpha;
% 
% parout.forces.AeroDyn_Forces = [A, N, T];
% parout.accelerations.body_acc = dV;
% 
% parout.air.rho = rho;
% parout.air.P = P;
% 
% parout.coeff.CD = CD;
% parout.coeff.XCP = XCP_value;





