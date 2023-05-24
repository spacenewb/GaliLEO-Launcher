
function [structural_mass_vect,eps_vect] = structural_masses (launcher, prop, data)

%% DATA EXTRAPOLATION FROM STRUCT 'DATA' :

mp1 = launcher.st1.mp;
mp2 = launcher.st2.mp;


m1 = launcher.st1.mp + launcher.st1.ms;
m2 = launcher.st2.mp + launcher.st2.ms;


c_star1 = launcher.st1.isp * data.orbit.g;
c_star2 = launcher.st2.isp * data.orbit.g;

p_c1 = prop.st1.Pc*1e-6; %missing 
p_c2 = prop.st2.Pc*1e-6; %missing 

t_burn1 = launcher.st1.tb;
t_burn2 = launcher.st2.tb;

eps_nozzle_1 = prop.st1.eps;
eps_nozzle_2 = prop.st2.eps; 

theta1 = prop.st1.alpha; %set
theta2 = prop.st2.alpha; %set

L1 = launcher.st1.Lcc;
L2 = launcher.st2.Lcc;
L3 = 0.4; % 1.95741 - 1.5; %[m] --> values from Hari's CAD model

rho_shell = data.str.rho_shell; %[kg/m^3]
rho_insulation = data.str.rho_ins; %[kg/m^3]
sigma = data.str.sigma;
sigma_hoop = data.str.sigma_hoop;
E = data.str.E;
FOS = data.str.FOS;

diam_ext = launcher.st1.C; %[m]

%% thickness with stress evaluation 

t_axial = (launcher.st1.T(1) + launcher.st1.m0*data.orbit.g)/(2*pi*(diam_ext/2) * (sigma/FOS));
Le = 2*(L1+L2+L3);
sigma_crit = (pi^2*E)/(Le/(diam_ext/2))^2;
t_buck = ((sigma_crit*FOS)/E) * (diam_ext/2) * (1/0.25); 
% t_buck = fsolve(fun,t_axial,data.optionsfsolve);
t_burst = (p_c1 * diam_ext) / (2*sigma_hoop);
 
shell_thick = sqrt( t_axial^2 + t_buck^2 + t_burst^2 );
insulation_thick = 1.2 * shell_thick; %[m]


% r = diam_ext/2;
% l = 1.5;
% fun_karman = @(x)(( r./sqrt(pi) ) .* sqrt( acos( 1- ( (2.*x)./l ) ) - ( sin(2.*acos( 1- ( (2.*x)./l ) )) ./ 2 ) )) .* sqrt(1 + ((1125899906842624*r.*(2./(l*(1 - ((2*x)/l - 1).^2).^(1/2)) - (2.*cos(2.*acos(1 - (2*x)./l)))/(l.*(1 - ((2.*x)./l - 1).^2).^(1/2))))./(3991211251234741.*(acos(1 - (2.*x)./l) - sin(2.*acos(1 - (2.*x)./l))./2).^(1/2))).^2);
% S_lat = integral(fun_karman,0,l) *2*pi;
% A_fairing = S_lat; %[m^2]

%% STAGE 1:

m_ext_shell_1 = (((diam_ext/2)^2)*pi*L1 - (((diam_ext-shell_thick)/2)^2)*pi*L1)*rho_shell + ((((diam_ext-shell_thick)/2)^2)*pi*L1 - (((diam_ext-shell_thick-insulation_thick)/2)^2)*pi*L1)*rho_insulation ;

m_nozzle_1 = 0.256*((10^(-4))*(((((mp1*c_star1)^1.2)*eps_nozzle_1^(0.3)))/((p_c1^0.8)*(t_burn1^0.6)*((tan(theta1))^0.4)))^0.917);

m_avionics_1 = 10*((m1)^0.361);

m_wiring_1 = 0; %1.508*(L1^0.25)*sqrt(m1);

m_recovery_1 = 100; %[kg]


m_struct_tot_stage1 = m_ext_shell_1 + m_nozzle_1 + m_avionics_1 + m_wiring_1 +m_recovery_1;
ms1 = m_struct_tot_stage1;

eps1 = ms1/(mp1+ms1);



%% STAGE 2:

m_ext_shell_2 = (((diam_ext/2)^2)*pi*L2 - (((diam_ext-shell_thick)/2)^2)*pi*L2)*rho_shell + ((((diam_ext-shell_thick)/2)^2)*pi*L2 - (((diam_ext-shell_thick-insulation_thick)/2)^2)*pi*L2)*rho_insulation ;

m_nozzle_2 = 0.256*((10^(-4))*(((((mp2*c_star2)^1.2)*eps_nozzle_2^(0.3)))/((p_c2^0.8)*(t_burn2^0.6)*((tan(theta2))^0.4)))^0.917);

m_avionics_2 = 10*((m2)^0.361);

m_wiring_2 = 0; % 1.508*(L2^0.25)*sqrt(m2);

m_struct_tot_stage2 = m_ext_shell_2 + m_nozzle_2 + m_avionics_2 + m_wiring_2 ;
ms2 = m_struct_tot_stage2;

eps2 = ms2/(mp2+ms2);



% %% STAGE 3:
% 
% m_ext_shell_3 = (((diam_ext/2)^2)*pi*L3 - (((diam_ext-shell_thick)/2)^2)*pi*L3)*rho_shell + ((((diam_ext-shell_thick)/2)^2)*pi*L3 - (((diam_ext-shell_thick-insulation_thick)/2)^2)*pi*L3)*rho_insulation ;
% 
% m_nozzle_3 = 0.256*((10^(-4))*(((((mp3*c_star3)^1.2)*eps_nozzle_3^(0.3)))/((p_c3^0.8)*(t_burn3^0.6)*((tan(theta3))^0.4)))^0.917);
% 
% m_fairing = 4.95*((A_fairing)^1.15) * 0.5;
% 
% 
% m_tanks = prop.st3.M_tank + prop.st3.Mpg_tank ; % missing 
% 
% m_avionics_3 = 10*((m3)^0.361);
% 
% m_wiring_3 = 1.508*(L3^0.25)*sqrt(m3);
% 
% m_T_struct_3 = 2.55*(10^(-4))*T3;
% 
% m_struct_tot_stage3 = m_ext_shell_3 +m_fairing+ m_nozzle_3 + m_avionics_3 + m_wiring_3 + m_T_struct_3 + m_tanks;
% ms3 = m_struct_tot_stage3;
% 
% eps3 = ms3/(mp3+ms3);

%% FINAL VECTOR:

structural_mass_vect = [ms1 ms2];
eps_vect = [eps1 eps2];

end