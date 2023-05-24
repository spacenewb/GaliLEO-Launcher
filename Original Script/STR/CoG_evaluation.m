function [x_cg,launcher] = CoG_evaluation(launcher, data)
% I consider a cylinder for each stage and a cone for the fairing
% first I evaluate the CoG considering each stage both full and empty and interpolating linearly
% then I evaluate the CoG for each second of burning considering a mean mass flow rate for each second
% for the third stage I consider half of the structural mass in the cylinder and half in the cone

g0 = data.orbit.g;

m1 = launcher.st1.mp + launcher.st1.ms;
m2 = launcher.st2.mp + launcher.st2.ms;
m3 = launcher.st3.mp + launcher.st3.ms;
m3_s = launcher.st3.ms;
mpay = launcher.st3.pay;

TW3_ratio = launcher.st3.T(1) / ((m3-m3_s/2)*g0);

Isp1 = launcher.st1.isp;
Isp2 = launcher.st2.isp;
Isp3 = launcher.st3.isp;

tb1 = round(launcher.st1.tb);
tb2 = round(launcher.st2.tb);
tb3 = round(launcher.st3.tb);

L1 = launcher.st1.L;
L2 = launcher.st2.L;
L_fairing = launcher.st3.fairing;
L3 = launcher.st3.L - L_fairing;
Ltot = L1 + L2 + L3 + L_fairing;

% Let's consider for smplicity the fairing as a cone with mass equals to half of
% the structural mass of third stage and CoG at 3/4*L from the tip of the nose

x_cg1_stat = Ltot - L1/2;       %from the tip of the nose
x_cg2_stat = Ltot - L1 - L2/2;
x_cg3_stat = (((L_fairing+L3/2)*(m3-m3_s/2))+(3/4*L_fairing*(mpay+m3_s/2)))/(m3+mpay);
x_cg3_stat_engine = Ltot - L1 - L2 - L3/2;
x_cg3_stat_payload = Ltot - L1 - L2 - L3 - L_fairing/4;
 
x_cg_tot(1) = Ltot - (m1*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + (m3+mpay)*(Ltot - x_cg3_stat))./(m1 + m2 + m3 + mpay);
x_cg_tot(2) = Ltot - (launcher.st1.ms*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + (m3+mpay)*(Ltot - x_cg3_stat))./(launcher.st1.ms + m2 + m3 + mpay);
x_cg_tot(3) = Ltot - (m2*(Ltot-x_cg2_stat) + (m3+mpay)*(Ltot-x_cg3_stat))/(m2 + m3 + mpay);
x_cg_tot(4) = Ltot - (launcher.st2.ms*(Ltot-x_cg2_stat) + (m3+mpay)*(Ltot-x_cg3_stat))/(launcher.st2.ms + m3 + mpay);
x_cg_tot(5) = Ltot - ((m3+mpay)*(Ltot-x_cg3_stat))/(m3 + mpay);
x_cg_tot(6) = Ltot - (m3_s/2*(Ltot-x_cg3_stat_engine) + (m3_s/2+mpay)*(Ltot-x_cg3_stat_payload))/(m3_s + mpay);

tb_vect = [0,tb1,tb1,tb1+tb2,tb1+tb2,tb1+tb2+tb3];

% figure
% plot(tb_vect,x_cg_tot, '-o')
% grid on

x_cg.x_cg_st1 = x_cg1_stat;
x_cg.x_cg_st2 = x_cg2_stat;
x_cg.x_cg_st3 = x_cg3_stat;
x_cg.x_cg_in = x_cg_tot(1);
x_cg.x_cg_vect_linear = x_cg_tot;

launcher.st1.xcgf = x_cg_tot(1);
launcher.st1.xcge = x_cg_tot(2);
launcher.st2.xcgf = x_cg_tot(3);
launcher.st2.xcge = x_cg_tot(4);
launcher.st3.xcgf = x_cg_tot(5);
launcher.st3.xcge = x_cg_tot(6);

%% prova non lineare

T1_vect = linspace(launcher.st1.T(1),launcher.st1.T(2),tb1);
mdot1_p = T1_vect./(Isp1*g0);
m1_vect = zeros(1,tb1+1);
m1_vect(1) = m1;
for i=2:tb1+1
    m1_vect(i) = m1_vect(i-1)-mdot1_p(i-1);
end
T2_vect = linspace(launcher.st2.T(1),launcher.st2.T(2),tb2);
mdot2_p = T2_vect./(Isp2*g0);

% no time accounted for the separation of the stages
m2_vect = zeros(1,tb2+1);
m2_vect(1) = m2;
for i=2:tb2+1
    m2_vect(i) = m2_vect(i-1)-mdot2_p(i-1);
end
m3_vect = (m3-m3_s/2).*(1-(TW3_ratio).*(linspace(0,tb3,tb3+1)./Isp3));

x_cg1_vect = Ltot - (m1_vect.*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + (m3+mpay)*(Ltot - x_cg3_stat))./(m1_vect + m2 + m3 + mpay);
x_cg2_vect = Ltot - (m2_vect.*(Ltot-x_cg2_stat) + (m3+mpay)*(Ltot-x_cg3_stat))./(m2_vect + m3 + mpay);
x_cg3_vect = Ltot - (m3_vect.*(Ltot-x_cg3_stat_engine) + (m3_s/2+mpay)*(Ltot-x_cg3_stat_payload))./(m3_vect + m3_s/2 + mpay);

tb1_vect = linspace(0, tb1, tb1+1);
tb2_vect = linspace(tb1, tb1+tb2, tb2+1);
tb3_vect = linspace(tb1+tb2, tb1+tb2+tb3, tb3+1);

% hold on
% plot(tb1_vect, x_cg1_vect,'-o', tb2_vect, x_cg2_vect,'-o', tb3_vect(1:3:end), x_cg3_vect(1:3:end), '-o')
% legend('linear','using mdot')

x_cg.x_cg_vect = [x_cg1_vect, x_cg2_vect, x_cg3_vect];

end