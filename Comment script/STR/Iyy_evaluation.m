function [launcher] = Iyy_evaluation(launcher, data)
% for semplicity let's consider the entire rocket as a cylinder, both for the stages either for the fairing
% first evaluation of the Iyy considering each stage both full and empty and interpolating linearly
% then evaluatio of the Iyy for each second of burning considering a mean mass flow rate for each second

g0 = data.orbit.g;  %gravity acceleration [m/s^2] 
stage = fieldnames(launcher);

% mass [kg]
m1 = launcher.st1.mp + launcher.st1.ms;
m2 = launcher.st2.mp + launcher.st2.ms;
m3 = launcher.st3.mp + launcher.st3.ms;
m3_s = launcher.st3.ms;     % structural mass stage 3 [kg]
mpay = launcher.st3.pay;    % payload mass [kg]

TW3_ratio = launcher.st3.T(1) / ((m3-m3_s/2)*g0);   % thrust over weight ratio

% Spacific Impulse [s]
Isp1 = launcher.st1.isp;   
Isp2 = launcher.st2.isp;  
Isp3 = launcher.st3.isp;   

% vector of burning time [s]
tb1 = round(launcher.st1.tb);   
tb2 = round(launcher.st2.tb);
tb3 = round(launcher.st3.tb);

% Thrust profile is regressive
% Let's consider a vector of Thrust linearly interpolated and compute the
% mass flow rate for each second
% Then compute the mass of the stage subtracting the
% relative mass flow rate for each second

T1_vect = linspace(launcher.st1.T(1),launcher.st1.T(2),tb1);    %Thrust vector stage 1 [N]
mdot1_p = T1_vect./(Isp1*g0);   %mass flow rate vector [kg/s]

% mass vector stage 1 [kg]
m1_vect = zeros(1,tb1+1);
m1_vect(1) = m1;
for i=2:tb1+1
    m1_vect(i) = m1_vect(i-1)-mdot1_p(i-1);
end

T2_vect = linspace(launcher.st2.T(1),launcher.st2.T(2),tb2);    %Thrust vector stage 2 [N]
mdot2_p = T2_vect./(Isp2*g0);   %mass flow rate vector [kg/s]

% no time accounted for the separation of the stages

% mass vector stage 2 [kg]
m2_vect = zeros(1,tb2+1);
m2_vect(1) = m2;
for i=2:tb2+1
    m2_vect(i) = m2_vect(i-1)-mdot2_p(i-1);
end

% Thrust is constant for stage 3
% mass vector stage 3 [kg]
m3_vect = (m3-m3_s/2).*(1-(TW3_ratio).*(linspace(0,tb3,tb3+1)./Isp3));

% mass vector for each stage [kg]
launcher.st1.m_vect = m1_vect + m2 + m3 + mpay;
launcher.st2.m_vect = m2_vect + m3 + mpay;
launcher.st3.m_vect = m3_vect + m3_s/2 + mpay;

% Computation of moment of inertia for each second [kg*m^2]
Iyy_linear = [];
for i = 1:3
    launcher.(stage{i}).Iyy_vect = (1/12) * launcher.(stage{i}).m_vect * (3*(launcher.(stage{i}).C/2)^2 + launcher.(stage{i}).L_stage^2);
    % vector of moment of inertia for linear interpolation [kg*m^2]
    Iyy_linear = [Iyy_linear,launcher.(stage{i}).Iyyf,launcher.(stage{i}).Iyye];
end

% time vectors for each stage [s]
tb1_vect = linspace(0, tb1, tb1+1);
tb2_vect = linspace(tb1, tb1+tb2, tb2+1);
tb3_vect = linspace(tb1+tb2, tb1+tb2+tb3, tb3+1);

% vector of each stage separation step [s]
tb_vect = [0,tb1,tb1,tb1+tb2,tb1+tb2,tb1+tb2+tb3];

figure
plot(tb_vect,Iyy_linear, '-x')
grid on
hold on
plot(tb1_vect,launcher.st1.Iyy_vect,'-o', tb2_vect,launcher.st2.Iyy_vect,'-o', tb3_vect,launcher.st3.Iyy_vect, '-o')
legend('linear','using mdot')
title('Moment of Inertia Iyy')
end