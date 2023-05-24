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
Is3 = 310.6; %[s]

%Exhaust velocities:
c1 = Is1*g/1000; %[km/s]
c2 = Is2*g/1000; %[km/s]
c3 = Is3*g/1000; %[km/s]

%structural coefficients:
eps1 = 0.09;
eps2 = 0.10;
eps3 = 0.20;

%% TWO STAGE OPTIMIZATION:

c = [c1 c2]'; %effective exhaust velocities vector
eps = [eps1 eps2]'; %structural coefficients vector

% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DeltaV_tot - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve

lambda0 = 1;
lambda = fsolve(fun, lambda0); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0
lambda_two = lambda;
% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage
mr_two = mr;

m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m_pay); %step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m2+m_pay); %step mass of first stage

mi_tot = m1+m2+m_pay;
mi_tot_two = m1+m2+m_pay;

m = [m1 m2];
m_two = [m1 m2];

% STRUCTURAL MASSES:

ms2 = m2*eps2;
ms1 = m1*eps1;

ms = [ms1 ms2]';
ms_two = [ms1 ms2]';

% PROPELLANT MASSES:

mp2 = m2 - ms2;
mp1 = m1 - ms1;

mp = [mp1 mp2]';
mp_two = [mp1 mp2]';

% PLOTS:

n_stages = [1:1:2];

figure(1)
plot (n_stages, m, 'g-');
hold on
grid on
plot (n_stages, ms, 'b-');
plot (n_stages, mp, 'r-');
title ('Two stages masses')
xlabel('Number of stages')
ylabel('Mass [kg]')
hold off

% DeltaVs FOR EACH STAGE:

DeltaV1 = c1*log(mr(1));
DeltaV2 = c2*log(mr(2));

DeltaV_stages = [DeltaV1 DeltaV2];
DeltaV_stages_two = DeltaV_stages;

%% THREE STAGE OPTIMIZATION:

c = [c1 c2 c3]'; %effective exhaust velocities vector
eps = [eps1 eps2 eps3]'; %structural coefficients vector

% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DeltaV_tot - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)- c3*log(lambda*c3-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve

lambda0 = 1;
lambda = fsolve(fun, lambda0); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0
lambda_three = lambda;

% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage
mr_three = mr;

m3 = ((mr(3)-1)/(1-eps3.*mr(3)))*m_pay; %step mass of third stage
m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m3+m_pay); %step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m3+m2+m_pay); %step mass of first stage

mi_tot = m1+m2+m3+m_pay;
mi_tot_three = mi_tot;

m = [m1 m2 m3];
m_three = m;

% STRUCTURAL MASSES:

ms3 = m3*eps3;
ms2 = m2*eps2;
ms1 = m1*eps1;

ms = [ms1 ms2 ms3]';
ms_three = ms;

% PROPELLANT MASSES:

mp3 = m3 - ms3;
mp2 = m2 - ms2;
mp1 = m1 - ms1;

mp = [mp1 mp2 mp3]';
mp_three = mp;

% PLOTS:

n_stages3 = [1:1:3];

figure(2)
plot (n_stages3, m, 'g-');
hold on
grid on
plot (n_stages3, ms, 'b-');
plot (n_stages3, mp, 'r-');
title ('Three stages masses')
xlabel('Number of stages')
ylabel('Mass [kg]')
hold off

% DeltaVs FOR EACH STAGE:

DeltaV1 = c1*log(mr(1));
DeltaV2 = c2*log(mr(2));
DeltaV3 = c3*log(mr(3));

DeltaV_stages = [DeltaV1 DeltaV2 DeltaV3];
DeltaV_stages_three = DeltaV_stages;

%% CONFRONTATION:
n_stages3 = [3 2 1];
n_stages2 = [2 1];

figure(3)
plot (n_stages3, DeltaV_stages_three, 'b-');
hold on
grid on
plot (n_stages2, DeltaV_stages_two, 'r-');
title ('Stages DeltaV confrontation')
xlabel('Number of stages')
ylabel('DeltaV [km/s]')
hold off


figure(4)
plot (n_stages3, m_three, 'b-');
hold on
grid on
plot (n_stages2, m_two, 'r-');
title ('Step masses confrontation')
xlabel('Number of stages')
ylabel('Step mass [kg]')
hold off

figure(5)
plot (n_stages3, ms_three, 'b-');
hold on
grid on
plot (n_stages2, ms_two, 'r-');
title ('Structural masses confrontation')
xlabel('Number of stages')
ylabel('Structural mass [kg]')
hold off

figure(6)
plot (n_stages3, mp_three, 'b-');
hold on
grid on
plot (n_stages2, mp_two, 'r-');
title ('Propellant masses confrontation')
xlabel('Number of stages')
ylabel('Propellant mass [kg]')
hold off

%% THRUST TO WEIGTH RATIOS:
k = [0.7:0.1:1.3];
s = length(k);

for j=1:s
    
T1 = k(j)*mi_tot*9.81;
tw1_start = T1/(mi_tot*9.81);
tw1_end = T1/((mi_tot-mp1)*9.81);

T2 = k(j)*(m2+m3+m_pay)*9.81;
tw2_start = T2/(m2*9.81);
tw2_end = T2/((m2+m3+m_pay-mp2)*9.81);

T3 = k(j)*(m3+m_pay)*9.81;
tw3_start = T3/(m3*9.81);
tw3_end = T3/((m3+m_pay-mp3)*9.81);

tw1 = [tw1_start tw1_end];
tw2 = [tw2_start tw2_end];
tw3 = [tw3_start tw3_end];
int = [1:2];

tw = [tw1_start tw1_end tw2_start tw2_end tw3_start tw3_end];
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
plot(int,tw3,'g-')
legend('First stage','Second stage','Third stage')
hold off
end
