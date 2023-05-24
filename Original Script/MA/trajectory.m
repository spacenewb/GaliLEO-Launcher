function [Y_tot,T_tot] = trajectory(launcher, site, N, orbit,options, sol)

x0 = 0;
z0 = site.z0;
V0 = 1;
gamma0 = deg2rad(site.OMEGA);

launcher.st3.tb = sol(1);

chi = sol(2:end);
stage = fieldnames(launcher);
k = 0;
t0_stage = 0;
Y_tot = [];
T_tot = [];
acc = [];

for i = 1:numel(stage)
    
    m0 = launcher.(stage{i}).m0;
    Y0 = [x0 z0 V0 gamma0 m0];
    tf = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf];
    
    dt = (tf-t0_stage)/N;

for j = 1:N 
   
    t0 = (j-1)*dt + t0_stage;
    tf = j*dt + t0_stage;
    [T, Y] = ode45(@ascent_opt, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi(j+k));
    Y0 = Y(end,:);
    Y_tot = [Y_tot; Y];
    T_tot = [T_tot; T];
    [parout] = RecallOde(@ascent_opt, T, Y, launcher.(stage{i}), site, orbit, chi(j+k));
    acc = [acc, parout.accelerations.body_acc];

end
    k = j;
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    t0_stage = tf;
 
   
end


%%

ind1 = find(T_tot==launcher.st1.tb);
ind2 = find(T_tot==(launcher.st1.tb + launcher.st2.tb));

figure
plot(T_tot,acc/9.81,'b','Linewidth',2)
title('Acceleration vs time')
xlabel('time [s]');
ylabel('a [g]');
grid on
hold on
scatter(T_tot(ind1(1)),acc(ind1(1))/9.81,60,'LineWidth',1.5);
scatter(T_tot(ind2(1)),acc(ind2(1))/9.81,60,'LineWidth',1.5);
legend('a','Staging 1','Staging 2')

