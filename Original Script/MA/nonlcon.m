function [c, ceq] = nonlcon(x, launcher, site, N,target,orbit,options)

x0 = 0;
z0 = site.z0;
V0 = 1;
gamma0 = deg2rad(site.OMEGA);

launcher.st3.tb = x(1);

chi = x(2:end);
stage = fieldnames(launcher);
k = 0;
t0_stage = 0;
Y_tot = [];


for i = 1:numel(stage)
    
    m0 = launcher.(stage{i}).m0;
    Y0 = [x0 z0 V0 gamma0 m0];
    tf = launcher.(stage{i}).tb + t0_stage;
    launcher.(stage{i}).T_time = [t0_stage tf];
    
    dt = (tf-t0_stage)/N;
    
for j = 1:N 
   
    t0 = (j-1)*dt + t0_stage;
    tf = j*dt + t0_stage;
    [~, Y] = ode45(@ascent_opt, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi(j+k));
    Y0 = Y(end,:);
    Y_tot = [Y_tot;Y];

end

    k = j;
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    t0_stage = tf;
    

end


% Event constraint
[~,~,~,rho] =  atmoscoesa(site.z0 + Y_tot(:,2),'none');
q_dyn = 0.5 .* ((Y_tot(:,3).^2) .* rho);
z = Y(end,2);
V = Y(end,3);
vx = V*cos(Y(end,4));
vz = V*sin(Y(end,4));


ceq = [target.vx - vx, target.vz - vz, target.z - z];
c = max(q_dyn) - 50e3;

