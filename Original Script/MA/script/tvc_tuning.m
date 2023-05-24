%%% TVC TUNING %%%
function tvc_tuning(launcher)

%% TVC 1
T = launcher.st1.T;
L = launcher.st1.L_stage;
x_cg = [launcher.st1.xcgf launcher.st1.xcge];
Iyy = [launcher.st1.Iyyf, launcher.st1.Iyye];

% Iyy = (Iyy(1)+Iyy(2))/2;
% T = (T(1)+T(2))/2;
% x_cg = (x_cg(1)+x_cg(2))/2;
T = T(1);
x_cg = x_cg(1);
Iyy = Iyy(1);

s = tf('s');
G = T*(x_cg-L)/Iyy/s^2;

% PD

ts = 0.05;
wc = 4.6/ts;
Ti = 20/wc;
Td = Ti/4;
LL = -(1+Td*s+1/(s*Ti))*G;
%sisotool(LL)
%pidtune(G,'PD')
%pidTuner(G,'PID')

% %% TVC 3
% 
% T = launcher.st3.T;
% L = launcher.st3.L_stage;
% x_cg = [launcher.st3.xcgf launcher.st3.xcge];
% Iyy = [launcher.st3.Iyyf, launcher.st3.Iyye];
% 
% Iyy = (Iyy(1)+Iyy(2))/2;
% T = (T(1)+T(2))/2;
% x_cg = (x_cg(1)+x_cg(2))/2;
% 
% s = tf('s');
% G = T*(x_cg-L)/Iyy/s^2;
% 
% % PD
% 
% ts = 0.05;
% wc = 4.6/ts;
% Ti = 20/wc;
% Td = Ti/4;
% LL = -(1+Td*s+1/(s*Ti))*G;
% %sisotool(LL)
% %pidtune(G,'PID')
% pidTuner(G,'PID')

