close all
clear all
clc

%% DATA:

g = 9.81; %[m/s] GRAVITY ACCELERATION ON EARTH
m_pay = 50; %[kg] PAYLOAD MASS
DeltaV_tot = 8.4586; %[km/s] TOTAL DELTA V NEEDED

%Specific impulses:
Is1 = 250.3; %[s]
Is2 = 266.3; %[s]
Is3 = 266.3; %[s]
Is4 = 310.6; %[s]

%Exhaust velocities:
c1 = Is1*g/1000; %[km/s]
c2 = Is2*g/1000; %[km/s]
c3 = Is3*g/1000; %[km/s]
c4 = Is4*g/1000; %[km/s]

c = [c1 c2 c3 c4]'; %effective exhaust velocities vector

%structural coefficients:
eps1 = 0.10;
eps2 = 0.15;
eps3 = 0.20;
eps4 = 0.25;

eps = [eps1 eps2 eps3 eps4]'; %structural coefficients vector

%% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DeltaV_tot - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)- c3*log(lambda*c3-1)- c4*log(lambda*c4-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve

lambda0 = 1;
lambda = fsolve(fun, lambda0); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0

%% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage

m4 = ((mr(4)-1)/(1-eps4.*mr(4)))*m_pay; %step mass of fourth stage
m3 = ((mr(3)-1)/(1-eps3.*mr(3)))*(m4+m_pay); %step mass of third stage
m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m4+m3+m_pay); %step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m4+m3+m2+m_pay); %step mass of first stage

mi_tot = m1+m2+m3+m4+m_pay;
m = [m1 m2 m3 m4];

%% STRUCTURAL MASSES:

ms4 = m4*eps4;
ms3 = m3*eps3;
ms2 = m2*eps2;
ms1 = m1*eps1;

ms = [ms1 ms2 ms3 ms4]';
%% PROPELLANT MASSES:

mp4 = m4 - ms4;
mp3 = m3 - ms3;
mp2 = m2 - ms2;
mp1 = m1 - ms1;

mp = [mp1 mp2 mp3 mp4]';

%% PLOTS:

n_stages = [1 2 3 4];

figure(1)
plot (n_stages, m, 'g-');
hold on
grid on
plot (n_stages, ms, 'b-');
plot (n_stages, mp, 'r-');
legend ('total masses','structural masses','propellant masses');
hold off

%% DeltaVs FOR EACH STAGE:

DeltaV1 = c1*log(mr(1));
DeltaV2 = c2*log(mr(2));
DeltaV3 = c3*log(mr(3));
DeltaV4 = c4*log(mr(4));

DeltaV_stages = [DeltaV1 DeltaV2 DeltaV3 DeltaV4];


