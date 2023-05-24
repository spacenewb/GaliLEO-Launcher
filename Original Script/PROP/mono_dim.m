

Tmax = 2000; % N valore sparato a caso
Pc = 10*10^5; % Pa
MM = 22.256;
%k = 1.1834;
T0 = 2131.11; % [K]
Ru = 8314.46;
Cp = 2.2125;
Cv = Cp - Ru/(MM*1000);
k = Cp/Cv; %better value for k
g0 = 9.81;
rho = 1357; % [kg/m^3]
mprop = 56.3048; % [kg]
alpha = 15; % [deg]
beta = 45; % [deg]
alpha = deg2rad(alpha);
beta = deg2rad(beta);
M_c = 0.2;

[stage_3] = third_stage(Tmax, Pc, MM, k, rho, mprop);
