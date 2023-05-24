function plot_sim(traj,data)        
Y = traj.Y;
T = traj.T;
% T1a = traj.T1a;
% Y1a = traj.Y1a;
% T1d = traj.T1d;
% Y1d = traj.Y1d;

figure
plot(T,Y(:,2)/1000,'b','Linewidth',2)
title('Altitude Profile vs Time')
xlabel('time [s]');
ylabel('z [km]');
grid on
hold on
% plot(T1a,Y1a(:,2)/1000,'r','Linewidth',2)
% plot(T1d,Y1d(:,3)/1000,'r','Linewidth',2)


figure
plot(Y(:,1)/1000,Y(:,2)/1000,'b','Linewidth',2)
title('Trajectory')
xlabel('x [km]');
ylabel('z [km]');
grid on
axis equal
hold on
% plot(Y1a(:,1)/1000,Y1a(:,2)/1000,'r','Linewidth',2)



figure
plot(T,Y(:,3),'b','Linewidth',2)
title('|V|')
xlabel('t [s]');
ylabel('V [m/s]');
grid on
hold on
% plot(T1a,Y1a(:,3),'r','Linewidth',2)


figure
plot(T,Y(:,3).*cos(Y(:,4)),'r','Linewidth',2)
title('Inertial velocities')
xlabel('t [s]');
ylabel('[m/s]');
grid on
hold on 
plot(T,Y(:,3).*sin(Y(:,4)),'b','Linewidth',2)
legend('vx','vz')


figure
plot(T,rad2deg(Y(:,4)),'b','Linewidth',2)
title('Attitude')
xlabel('time [s]');
ylabel('theta [°]');
grid on

% Mach 
figure
plot(T,traj.interp.M,'b','Linewidth',2)
title('Mach')
xlabel('time [s]');
ylabel('M');
grid on
hold on 
% plot(T1a,traj.bal.interp.M,'r','Linewidth',2)
% legend('Path to orbit','First stage')

% Stagnation temperature
index = Y(:,2) < 60000;
[Tamb, ~, ~, ~] = atmoscoesa(Y(index,2),'none');
Ttot = Tamb'.*(1+0.2.*traj.interp.M(index).^2);


figure
plot(T(index),Ttot,'b','Linewidth',2)
title('Stagnation temperature')
xlabel('time [s]');
ylabel('T_{tot} [k]');
grid on




% dynamic pressure
q_dyn = 0.5 .* ((Y(:,3).^2)' .* traj.air.rho);
% q_dyn_a = 0.5 .* ((Y1a(:,3).^2)' .* traj.bal.air.rho);
figure
plot(T,q_dyn,'b','Linewidth',2)
title('Dynamic pressure')
xlabel('time [s]');
ylabel('Dynamic pressure [Pa]');
grid on
hold on
% plot(T1a,q_dyn_a,'r','Linewidth',2)

figure
plot(T,traj.accelerations/9.81,'b','Linewidth',2)
title('x-body acc')
xlabel('time [s]');
ylabel('x.acc [g]');
grid on

figure
plot(T,traj.aero_force(:,1),'b','Linewidth',2)
title('Forces')
xlabel('time [s]');
ylabel('[N]');
grid on
hold on 
plot(T,traj.aero_force(:,3),'r','Linewidth',2)
legend('Drag','Thrust')


figure
xnorth = cosd(data.site.azimuth).*Y(:,1);
yeast = sind(data.site.azimuth).*Y(:,1);
% index = (T < T1a(1));
% xnorth_st1 = [cosd(data.site.azimuth).*Y(index,1); cosd(data.site.azimuth).*Y1a(:,1); Y1d(:,1)];
% yeast_st1 = [sind(data.site.azimuth).*Y(index,1); sind(data.site.azimuth).*Y1a(:,1); Y1d(:,2)];
[lat,lon,h] = ned2geodetic(xnorth, yeast,-Y(:,2), data.site.lat0, data.site.lon0, -data.site.z0, wgs84Ellipsoid);
% [lat_st1,lon_st1] = ned2geodetic(xnorth_st1, yeast_st1, 0, data.site.lat0, data.site.lon0, 0, wgs84Ellipsoid);
geoplot(lat,lon,'r','Linewidth', 2)
% hold on 
% geoplot(lat_st1,lon_st1,'b','Linewidth', 2)
geobasemap streets
legend('Path to orbit','First stage path')

% uif = uifigure;
% g = geoglobe(uif);
% geoplot3(g,lat,lon,h,'r','Linewidth', 2)

figure
plot(T,traj.coeff.CD,'b','Linewidth',2)
title('CD')
xlabel('time [s]');
ylabel('CD');
grid on

figure
plot(T,rad2deg(traj.interp.alpha),'b','Linewidth',2)
title('AoA')
xlabel('time [s]');
ylabel('alpha');
grid on