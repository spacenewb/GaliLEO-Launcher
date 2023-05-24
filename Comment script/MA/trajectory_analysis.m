function trajectory_analysis(launcher,orbit,site,data)

options_ascent = data.options_ascent;
options_descent = data.options_descent;
obj = load('orbit_path.mat');


%% Data
launcher.st0 = launcher.st1;
launcher.st0.m = launcher.st1.ms;
launcher.st0.aero.full = load('st1.mat');
launcher.st0.aero.empty = load('st1.mat');

%% Ascent First stage
% ballistic flight of first stage after burnout before drogue opening
launcher.st0.T_time = [0 launcher.st1.tb];
Y01 = obj.Y(obj.T==launcher.st1.tb,:)';
Y01 = Y01(:,1);
Y01(5) = launcher.st0.m;
Y01(6) = 0;
chi = 0;
[T1a, Y1a] = ode113(@rocket_dyn, [launcher.st1.tb 2000], Y01, options_ascent, launcher.st0, site, orbit, chi);
[traj.bal] = RecallOde(@rocket_dyn, T1a, Y1a, launcher.st0, site, orbit, chi);


%% Descent 
% descent with parachute

% initial state
t0 = T1a(end);
x0d = Y1a(end,1)*cosd(site.azimuth);
y0d = Y1a(end,1)*sind(site.azimuth);
z0d = Y1a(end,2);
vxd = Y1a(end,3)*cosd(site.azimuth)*cos(Y1a(end,4));
vyd = Y1a(end,3)*sind(site.azimuth)*cos(Y1a(end,4));
vzd = Y1a(end,3)*sin(Y1a(end,4));
Y01d = [x0d y0d z0d vxd vyd vzd];
T1d_tot = [];
Y1d_tot = [];

for i = 1 : length(data.para)
    

[T1d, Y1d] = ode45(@descent, [t0 6000], Y01d, options_descent, launcher.st0, site, orbit, data.para(i));
Y01d = Y1d(end,:);

t0 = T1d(end);
d_t(i) = T1d(end);

T1d_tot = [T1d_tot; T1d];
Y1d_tot = [Y1d_tot; Y1d];

end

ind1d = find(T1d_tot==d_t(1));
ind2d = find(T1d_tot==d_t(2));




% RECALL PLOT QUONTITIES
ind1 = find(obj.T==launcher.st1.tb);
ind2 = find(obj.T==(launcher.st1.tb + launcher.st2.tb));
T1 = obj.T(1:ind1(1));
T2 = obj.T(ind1(2):ind2(1));
T3 = obj.T(ind2(2):end);
Y1 = obj.Y(1:ind1(1),:);
Y2 = obj.Y(ind1(2):ind2(1),:);
Y3 = obj.Y(ind2(2):end,:);


launcher.st1.T_time = [0 T1(end)];
launcher.st2.T_time = [T1(end) T2(end)];
launcher.st3.T_time = [T2(end) T3(end)];

[traj.st1] = RecallOde(@ascent_opt, T1,Y1 , launcher.st1, site, orbit, chi);
[traj.st2] = RecallOde(@ascent_opt, T2, Y2, launcher.st2, site, orbit, chi);
[traj.st3] = RecallOde(@ascent_opt, T3, Y3, launcher.st3, site, orbit, chi);

traj.air.P = [traj.st1.air.P traj.st2.air.P traj.st3.air.P];
traj.air.rho = [traj.st1.air.rho traj.st2.air.rho traj.st3.air.rho];
traj.accelerations = [traj.st1.accelerations.body_acc, traj.st2.accelerations.body_acc, traj.st3.accelerations.body_acc];
traj.interp.M = [traj.st1.interp.M traj.st2.interp.M traj.st3.interp.M];
traj.interp.alpha = [traj.st1.interp.alpha traj.st2.interp.alpha traj.st3.interp.alpha];
traj.aero_force = [traj.st1.forces.AeroDyn_Forces; traj.st2.forces.AeroDyn_Forces; traj.st3.forces.AeroDyn_Forces];
traj.coeff.CD = [traj.st1.coeff.CD, traj.st2.coeff.CD, traj.st3.coeff.CD];




%% PLOT  

figure
plot(obj.T,obj.Y(:,2)/1000,'b--','Linewidth',2)
title('Altitude Profile vs Time')
xlabel('time [s]');
ylabel('Altitude [km]');
grid on
hold on
plot(T1a,Y1a(:,2)/1000,'m','Linewidth',2)
plot(T1d_tot,Y1d_tot(:,3)/1000,'r','Linewidth',2)
scatter(obj.T(ind1(1)),obj.Y(ind1(1),2)/1000,60,'LineWidth',1.5);
scatter(T1a(end),Y1a(end,2)/1000,60,'LineWidth',1.5);
scatter(T1d_tot(ind1d),Y1d_tot(ind1d,3)/1000,60,'LineWidth',1.5);
scatter(T1d_tot(ind2d),Y1d_tot(ind2d,3)/1000,60,'LineWidth',1.5);
scatter(T1d_tot(end),Y1d_tot(end,3)/1000,60,'LineWidth',1.5);
%plot(T2a,Y2a(:,2)/1000,'Linewidth',2)
legend('Ascent path','Stage 1 ballistic flight','Descent parachute','Staging 1','Opening 1° drogue','Opening 2° drogue','Opening main parachute','Landing')



figure
yyaxis left
plot(obj.Y(:,1)/1000,obj.Y(:,2)/1000,'b','Linewidth',2)
title('Trajectory and acceleration')
xlabel('Range [km]');
ylabel('Altitude [km]');
grid on
hold on
%plot(Y1a(:,1)/1000,Y1a(:,2)/1000,'r','Linewidth',2)
scatter(obj.Y(ind1(1),1)/1000,obj.Y(ind1(1),2)/1000,60,'LineWidth',1.5);
scatter(obj.Y(ind2(1),1)/1000,obj.Y(ind2(1),2)/1000,60,'LineWidth',1.5);
%plot(Y2a(:,1)/1000,Y2a(:,2)/1000,'Linewidth',2)
% plot(Y3a(:,1)/1000,Y3a(:,2)/1000,'Linewidth',2)
yyaxis right
plot(obj.Y(:,1)/1000,traj.accelerations/9.81,'r','Linewidth',2)
% title('Acceleration vs time')
% xlabel('time [s]');
ylabel('a [g]');
% grid on
% hold on
scatter(obj.Y(ind1(1),1)/1000,traj.accelerations(ind1(1))/9.81,60,'LineWidth',1.5);
scatter(obj.Y(ind2(1),1)/1000,traj.accelerations(ind2(1))/9.81,60,'LineWidth',1.5);
% legend('a','Staging 1','Staging 2')

legend('Trajectory','Staging 1','Staging 2','a','Staging 1','Staging 2')


figure
plot(obj.T,rad2deg(obj.Y(:,4)),'b','Linewidth',2)
title('Flight path angle vs time')
xlabel('t [s]');
ylabel('\gamma [°]');
grid on
hold on
scatter(obj.T(ind1(1)),rad2deg(obj.Y(ind1(1),4)),60,'LineWidth',1.5);
scatter(obj.T(ind2(1)),rad2deg(obj.Y(ind2(1),4)),60,'LineWidth',1.5);
legend('\gamma','Staging 1','Staging 2')



figure
plot(obj.T,obj.Y(:,3),'b','Linewidth',2)
title('Module of velocity vs time')
xlabel('t [s]');
ylabel('V [m/s]');
grid on
hold on
%plot(T1a,Y1a(:,3),'r','Linewidth',2)
scatter(obj.T(ind1(1)),obj.Y(ind1(1),3),60,'LineWidth',1.5);
scatter(obj.T(ind2(1)),obj.Y(ind2(1),3),60,'LineWidth',1.5);
legend('|V|','Staging 1','Staging 2')


% Mach 
figure
plot(obj.T,traj.interp.M,'b','Linewidth',2)
title('Mach')
xlabel('time [s]');
ylabel('M');
grid on
hold on 
% plot(T1a,traj.bal.interp.M,'r','Linewidth',2)
% legend('Path to orbit','First stage')


figure
plot(obj.T,traj.accelerations/9.81,'b','Linewidth',2)
title('Acceleration vs time')
xlabel('time [s]');
ylabel('a [g]');
grid on
hold on
scatter(obj.T(ind1(1)),traj.accelerations(ind1(1))/9.81,60,'LineWidth',1.5);
scatter(obj.T(ind2(1)),traj.accelerations(ind2(1))/9.81,60,'LineWidth',1.5);
legend('a','Staging 1','Staging 2')


figure
plot(obj.T,traj.aero_force(:,1),'b','Linewidth',2)
title('Forces')
xlabel('time [s]');
ylabel('[N]');
grid on
hold on 
plot(obj.T,traj.aero_force(:,3),'r','Linewidth',2)
%legend('Drag','Thrust')


% dynamic pressure
q_dyn = 0.5 .* ((obj.Y(:,3).^2)' .* traj.air.rho);
figure
yyaxis left
plot(obj.T,q_dyn,'b','Linewidth',2)
title('Dynamic pressure and velocity')
xlabel('Time [s]');
ylabel('Dynamic pressure [Pa]');
grid on
hold on
yyaxis right 
plot(obj.T,obj.Y(:,3),'r','Linewidth',2)
scatter(obj.T(ind1(1)),obj.Y(ind1(1),3),60,'LineWidth',1.5);
scatter(obj.T(ind2(1)),obj.Y(ind2(1),3),60,'LineWidth',1.5);
ylabel('Velocity [m/s]')
legend('q_{dyn}','|V|','Staging 1','Staging 2')
% plot(T1a,q_dyn_a,'r','Linewidth',2)

figure
grid on
hold on
yyaxis left
plot(obj.Y(:,1)/1000,obj.Y(:,3).*cos(obj.Y(:,4)),'b','Linewidth',2)
title('Velocity components')
xlabel('Range [km]');
ylabel('v_x [m/s]');
yyaxis right 
plot(obj.Y(:,1)/1000,obj.Y(:,3).*sin(obj.Y(:,4)),'r','Linewidth',2)
ylabel('v_z [m/s]')
% scatter(obj.T(ind1(1)),obj.Y(ind1(1),3),60,'LineWidth',1.5);
% scatter(obj.T(ind2(1)),obj.Y(ind2(1),3),60,'LineWidth',1.5);
legend('v_x','v_z')




figure();
xnorth = cosd(data.site.azimuth).*obj.Y(:,1);
yeast = sind(data.site.azimuth).*obj.Y(:,1);

index = (obj.T < T1a(1));
xnorth_st1_b = [cosd(data.site.azimuth).*obj.Y(index,1); cosd(data.site.azimuth).*Y1a(:,1); Y1d_tot(:,1)];
yeast_st1_b = [sind(data.site.azimuth).*obj.Y(index,1); sind(data.site.azimuth).*Y1a(:,1); Y1d_tot(:,2)];
[lat,lon,h] = ned2geodetic(xnorth, yeast,-obj.Y(:,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
[lat_st1,lon_st1] = ned2geodetic(xnorth(ind1(1)), yeast(ind1(1)),-obj.Y((ind1(1)),2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
% [lat_d1,lon_d1] = ned2geodetic(Y1a(end,1) .* cosd(data.site.azimuth),Y1a(end,1) .* sind(data.site.azimuth),-Y1a(end,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
% [lat_d2,lon_d2] = ned2geodetic(Y1d_tot(ind1d,1), Y1d_tot(ind1d,2),-Y1d_tot(ind1d,3), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
% [lat_d3,lon_d3] = ned2geodetic(Y1d_tot(ind2d,1), Y1d_tot(ind2d,2),-Y1d_tot(ind2d,3), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
[lat_dl,lon_dl] = ned2geodetic(Y1d_tot(end,1), Y1d_tot(end,2),-Y1d_tot(end,3), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
[lat_st2,lon_st2] = ned2geodetic(xnorth(ind2(1)), yeast(ind2(1)),-obj.Y((ind2(1)),2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
[lat_st3,lon_st3] = ned2geodetic(xnorth(end), yeast(end),-obj.Y(end,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
%[lat2,lon2,h2] = ned2geodetic(xnorth2, yeast2,Y2a(:,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
%[lat3,lon3,h3] = ned2geodetic(xnorth3, yeast3,Y3a(:,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
[lat_st1_b,lon_st1_b] = ned2geodetic(xnorth_st1_b, yeast_st1_b, 0, data.site.lat0, data.site.lon0, 0, wgs84Ellipsoid);
geoplot(lat,lon,'b','Linewidth', 2)
% geoplot(lat_st1_b,lon_st1_b,'b','Linewidth', 2)
hold on 
geoscatter(lat_st1,lon_st1,60,'LineWidth',1.5)
% geoscatter(lat_d1,lon_d1,60,'LineWidth',1.5)
% geoscatter(lat_d2,lon_d2,60,'LineWidth',1.5)
% geoscatter(lat_d3,lon_d3,60,'LineWidth',1.5)
% geoscatter(lat_dl,lon_dl,60,'LineWidth',1.5)
geoscatter(lat_st2,lon_st2,60,'LineWidth',1.5)
geoscatter(lat_st3,lon_st3,60,'LineWidth',1.5)
% geoplot(lat2,lon2,'Linewidth', 2)
% geoplot(lat3,lon3,'Linewidth', 2)
% R = 10e3 ;    % radius of circle 
% C = [Y1d_tot(end,1), Y1d_tot(end,2)] ;   % center of circle 
% th = linspace(0,2*pi,1000);
% x_circ  = C(1)+R*cos(th)'; 
% y_circ = C(2)+R*sin(th)';
% [lat_c,lon_c] = ned2geodetic(x_circ, y_circ ,0, data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
% geoplot(lat_c,lon_c,'r--','Linewidth',1.5)
geobasemap streets
legend('Ascent path','Staging 1','Staging 2','Perigee') %,'First stage path')
% legend('First stage path','Staging 1','Opening drogue 1','Opening drogue 2','Opening main','Landing','Safety circle') %,'First stage path')
hold off


x = obj.Y(:,1)/1000;
omega = x./orbit.Re;
z = obj.Y(:,2)/1000 + orbit.Re;

x_a = z.*cos(omega);
y_a = z.*sin(omega);
z_a =  zeros(length(x_a),1);

a = 6653;
e = 0.0188;
earthplot = figure();
title('Trajectory around Earth')
plot3(x_a,y_a,z_a,'Linewidth',2);
hold on
plotOrbit(6378+400,0,0,0,omega(end),pi,3*pi,0.01,orbit.mu,'r--','2D',false)
plotOrbit(a,e,0,0,omega(end),0,pi,0.01,orbit.mu,'k','2D',false)
scatter3(x_a(end),y_a(end),z_a(end),60,'Linewidth',1.5);
scatter3(-6467.15,-2029.34,0,60,'Linewidth',1.5);
plotPlanet(3, [0 0 0],earthplot);
xlabel('X [km]')
ylabel('Y [km]')
xlim([-8e3 8e3])
ylim([-8e3 8e3])

legend( 'Ascent path','Target orbit','Elliptical transfer orbit','First third stage cut-off','Final burning')










