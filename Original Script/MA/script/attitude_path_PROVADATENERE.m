clear all
close all
clc

load traj.mat
gamma = Y(:,4);
m = length(T);

%% fitting

f = fit(T,gamma,'smoothingspline');
gamma_obj = feval(f,T(1):0.1:T(end));

figure()
plot(T(1):0.1:T(end),gamma_obj,'r','LineWidth',3)
grid on
hold on
plot(T,gamma,'bo')


