clear all
clc
run MAIN

Y0(1) = 0;
Y0(2) = 0;
Y0(3) = 0;
Y0(4) = 0;
Y0(5) = deg2rad(90);
Y0(6) = 0;

tspan = [0:1:30];
wind = -30; %[m/s]
s = length(wind);
chi = deg2rad(0);
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8,'Events',@stop);


winds = wind;
[T,Y] = ode45(@start_dynamic, tspan, Y0, options, launcher, chi, winds);
    




