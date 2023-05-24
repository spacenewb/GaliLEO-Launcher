function [dY] = start_dynamic (t, Y, launcher,chi,wind)

x = Y(1);
z = Y(2);
vx = Y(3);
vz = Y(4);
theta = Y(5);
p = Y(6);

V = sqrt(vx^2+vz^2);
vxr = vx-wind;

if vxr > 10^(-6)
    alpha = theta-atan2(vz,vxr);
else
    alpha = 0;
end
  
%% ASSEMBLE DERIVATES
T = launcher.st1.T(1);
m = launcher.st1.m0;
Iyy = launcher.st1.Iyyf;
g = 9.81;
rho = 1.225;
Re = 6378;
r = z+Re;
CL = 3.3*alpha;
CMA = 18;
S = launcher.st1.S;
C = launcher.st1.C;
L_stage = launcher.st1.L_stage;
q_dyn = 0.5*rho*V^2;
N = q_dyn*S*CL;
M = q_dyn*S*C*CMA*alpha;
xcg = launcher.st1.xcgf;
Mtvc = T*sin(chi)*(xcg-L_stage);    
    
Fx = T*cos(theta+chi)-N*sin(theta);
Fy = T*sin(theta+chi)+N*cos(theta)-m*g;
dp = (M + Mtvc )/Iyy;

dY(1) = vx;
dY(2) = vz;
dY(3) = Fx/m;
dY(4) = Fy/m;
dY(5) = p;
dY(6) = dp;

dY = dY';

end