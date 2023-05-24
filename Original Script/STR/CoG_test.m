function [x_cg] = CoG_test(launcher, data)

[struct_mass_vect, ~] = structural_masses (launcher, data.prop, data);

ms=struct_mass_vect;
g0 = data.orbit.g;

m1 = launcher.st1.mp + launcher.st1.ms;
m2 = launcher.st2.mp + launcher.st2.ms;
m3 = launcher.st3.mp + launcher.st3.ms;

TW1_ratio = launcher.st1.tw_ratio;
TW2_ratio = launcher.st2.tw_ratio;
TW3_ratio = launcher.st3.T(1) / (m3*g0); %miss the payload mass?

Isp1 = launcher.st1.isp;
Isp2 = launcher.st2.isp;
Isp3 = launcher.st3.isp;

tb1 = launcher.st1.tb;
tb2 = launcher.st2.tb;
tb3 = launcher.st3.tb;

L1 = launcher.st1.Lcc; %is it the total lenght?
L2 = launcher.st2.Lcc;
L3 = 1.95741; %[m] --> values from Hari's CAD model
Ltot = L1 + L2 + L3;

x_cg1_stat = Ltot - L1/2;       %from the tip of the nose
x_cg2_stat = Ltot - L1 - L2/2;
x_cg3_stat = Ltot - L1 - L2 - L3/2;
x_cg_tot_0 = Ltot - (m1*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + m3*(Ltot - x_cg3_stat))./(m1 + m2 + m3);
x_cg_tot_1 = Ltot - (launcher.st1.ms*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + m3*(Ltot - x_cg3_stat))./(launcher.st1.ms + m2 + m3);
x_cg_tot_2 = Ltot - (m2*(Ltot-x_cg2_stat) + m3*(Ltot-x_cg3_stat))/(m2 + m3);
x_cg_tot_3 = Ltot - (launcher.st2.ms*(Ltot-x_cg2_stat) + m3*(Ltot-x_cg3_stat))/(launcher.st2.ms + m3);
x_cg_tot_4 = Ltot - (m3*(Ltot-x_cg3_stat))/m3;

tb_vect = [0,tb1,tb1,tb1+tb2,tb1+tb2,tb1+tb2+tb3];

x_cg_vect = [x_cg_tot_0,x_cg_tot_1,x_cg_tot_2,x_cg_tot_3,x_cg_tot_4];

figure
plot(tb_vect,[x_cg_vect,x_cg_tot_4], '-o')
grid on

x_cg.x_cg_st1 = x_cg1_stat;
x_cg.x_cg_st2 = x_cg2_stat;
x_cg.x_cg_st3 = x_cg3_stat;
x_cg.x_cg_tot = x_cg_tot_0;
x_cg.x_cg_vect = x_cg_vect;

%% prova non lineare

T1_vect = linspace(launcher.st1.T(1),launcher.st1.T(2),round(tb1));
mdot1_p = T1_vect./(Isp1*g0);
m1_vect = zeros(1,round(tb1)+1);
m1_vect(1) = m1;
for i=2:round(tb1+1)
    m1_vect(i) = m1_vect(i-1)-mdot1_p(i-1);
end
T2_vect = linspace(launcher.st2.T(1),launcher.st2.T(2),round(tb2));
mdot2_p = T2_vect./(Isp2*g0);
m2_vect = zeros(1,round(tb2)+1);
m2_vect(1) = m2;
for i=2:round(tb2+1)
    m2_vect(i) = m2_vect(i-1)-mdot2_p(i-1);
end
m3_vect = m3.*(1-(TW3_ratio).*(linspace(1,tb3,round(tb3))./Isp3));

x_cg1_vect = Ltot - (m1_vect.*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + m3*(Ltot - x_cg3_stat))./(m1_vect + m2 + m3);
x_cg2_vect = Ltot - (m2_vect.*(Ltot-x_cg2_stat) + m3*(Ltot-x_cg3_stat))./(m2_vect +m3);
x_cg3_vect = Ltot - (m3_vect.*(Ltot-x_cg3_stat))./m3_vect;

tb = round(tb1+1)+round(tb2+1)+round(tb3);
tb_vect = linspace(0,tb,tb);
tb_vect = tb_vect(1:3:end);
x_cg_vect = [x_cg1_vect,x_cg2_vect,x_cg3_vect];
x_cg_vect = x_cg_vect(1:3:end);

hold on
plot(tb_vect, x_cg_vect, '-o')

%% prova nuove masse strutturali

m1 = launcher.st1.mp + ms(1);
m2 = launcher.st2.mp + ms(2);
m3 = launcher.st3.mp + ms(3);
mpay = 50;

x_cg_tot_0 = Ltot - (m1*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + m3*(Ltot - x_cg3_stat))./(m1 + m2 + m3 );
x_cg_tot_1 = Ltot - (launcher.st1.ms*(Ltot - x_cg1_stat) + m2*(Ltot - x_cg2_stat) + m3*(Ltot - x_cg3_stat))./(launcher.st1.ms + m2 + m3);
x_cg_tot_2 = Ltot - (m2*(Ltot-x_cg2_stat) + m3*(Ltot-x_cg3_stat))/(m2 + m3);
x_cg_tot_3 = Ltot - (launcher.st2.ms*(Ltot-x_cg2_stat) + m3*(Ltot-x_cg3_stat))/(launcher.st2.ms + m3);
x_cg_tot_4 = Ltot - (m3*(Ltot-x_cg3_stat))/(m3);

tb_vect = [0,tb1,tb1,tb1+tb2,tb1+tb2,tb1+tb2+tb3];

x_cg_vect = [x_cg_tot_0,x_cg_tot_1,x_cg_tot_2,x_cg_tot_3,x_cg_tot_4,x_cg_tot_4];

hold on
plot(tb_vect,x_cg_vect, 'g-o')
legend('linear old masses','old masses','linear new masses')

end