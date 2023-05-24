function cost = dv_eval(x,launcher,orbit,site,N,options)

% cost = x(1);

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
T_tot = [];

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

end
    k = j;
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    t0_stage = tf;
 
   
end

   
    CD = 0.49;
    [~,~,~,rho] = atmoscoesa(Y_tot(:,2),'none');
    V = Y_tot(:,3);
    Re = orbit.Re * 1e3;
    drag = 0.5*rho.*V.^2.*launcher.st1.S.*CD;
    drag_loss = drag./Y_tot(:,5);
    index = not(isnan(drag_loss));
    dv_drag = trapz(T_tot(index),drag_loss(index));
    grav_loss = orbit.g .* ( Re./(Re + Y_tot(:,2)) ).^2 .* sin(Y_tot(:,4));
    dv_grav = trapz(T_tot,grav_loss);
    
    cost = dv_grav+dv_drag;




