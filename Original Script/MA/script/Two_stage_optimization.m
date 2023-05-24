close all
clear all
clc

%% DATA:

g = 9.81; %[m/s] GRAVITY ACCELERATION ON EARTH
m_pay = 50; %[kg] PAYLOAD MASS
DeltaV_tot = 10; %[km/s] TOTAL DELTA V NEEDED

%Specific impulses:
Is1 = 400; %[s]
Is2 = 350; %[s]

%Exhaust velocities:
c1 = Is1*g/1000; %[km/s]
c2 = Is2*g/1000; %[km/s]

c = [c1 c2]'; %effective exhaust velocities vector

%structural coefficients:
eps1 = 0.10;
eps2 = 0.15;

eps = [eps1 eps2]'; %structural coefficients vector

%% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DeltaV_tot - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve

lambda0 = 1;
lambda = fsolve(fun, lambda0); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0

%% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage

m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m_pay); %step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m2+m_pay); %step mass of first stage

mi_tot = m1+m2+m_pay;
m = [m1 m2];

%% STRUCTURAL MASSES:

ms2 = m2*eps2;
ms1 = m1*eps1;

ms = [ms1 ms2]';
%% PROPELLANT MASSES:

mp2 = m2 - ms2;
mp1 = m1 - ms1;

mp = [mp1 mp2]';

%% PLOTS:

n_stages = [1 2];

figure(1)
plot (n_stages, m, 'g-');
hold on
grid on
plot (n_stages, ms, 'b-');
plot (n_stages, mp, 'r-');
hold off

%% DeltaVs FOR EACH STAGE:

DeltaV1 = c1*log(mr(1));
DeltaV2 = c2*log(mr(2));

DeltaV_stages = [DeltaV1 DeltaV2];

%% THRUST TO WEIGTH RATIOS:
k = [0.7:0.1:1.1];
s = length(k);

for j=1:s
    
T1 = k(j)*mi_tot*9.81;
tw1_start = T1/(mi_tot*9.81);
tw1_end = T1/((mi_tot-mp1)*9.81);

T2 = k(j)*(m2+m_pay)*9.81;
tw2_start = T2/(m2*9.81);
tw2_end = T2/((m2+m_pay-mp2)*9.81);

tw1 = [tw1_start tw1_end];
tw2 = [tw2_start tw2_end];
int = [1:2];

tw = [tw1_start tw1_end tw2_start tw2_end];
n = length(tw);
tw_max = tw(1);
for i=2:n
    if (tw(n)>tw_max)
        tw_max = tw(n);
    end
end

nn = k(j);
figure (6+j)
plot(int,tw1,'b-')
title('Thrust to weigth ratio variation')
xlabel('start to end burning')
ylabel('T/W')
hold on 
grid on
plot(int,tw2,'r-')
legend('First stage','Second stage','Third stage')
hold off
end

