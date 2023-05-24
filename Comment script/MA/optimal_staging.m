function [launcher,DeltaV_stages] = optimal_staging(data)

orbit = data.orbit; 
launcher = data.launcher;

%% DATA:

% from the final orbit data contaneid in data.orbit the orbital parameters
% of the transfer orbit such perigee and apogee radius are here computed
% and with them the respective velocities:
rp = orbit.rp_e + orbit.Re;
ra = orbit.h + orbit.Re;
e = (ra-rp)/(ra+rp);
a = (ra+rp)/2;
p = a*(1-e^2);
vp = sqrt(orbit.mu/p)*(1+e); %transfer prbit perigee velocity
va = sqrt(orbit.mu/p)*(1-e); %transfer orbit apogee velocity

v_circ = sqrt( orbit.mu/(orbit.h+orbit.Re) ) ; %velocity of the final circular orbit
dv_orbit = v_circ - va; %DV for circularition impulse to move from the transfer orbit to the final one


dv_loss = orbit.losses/1000; %DV losses in [km/s]
DV = (vp + dv_orbit + dv_loss); %total DV as summation of the velocity to reach the perigee, the circularization DV and the losses
g = orbit.g; %[m/s] GRAVITY ACCELERATION ON EARTH

m_pay = launcher.st3.pay; %[kg] PAYLOAD MASS

%structural coefficients:
eps1 = launcher.st1.eps; %structural coefficent of stage 1
eps2 = launcher.st2.eps; %structural coefficent of stage 2
eps3 = launcher.st3.eps; %structural coefficent of stage 3

eps = [eps1 eps2 eps3]'; %structural coefficients vector

%Specific impulses:
Is1 = launcher.st1.isp; %[s] specific impulse of stage 1
Is2 = launcher.st2.isp; %[s] specific impulse of stage 2
Is3 = launcher.st3.isp; %[s] specific impulse of stage 3

%Exhaust velocities:
c1 = Is1*g/1000; %[km/s] exhaust velocity of stage 1
c2 = Is2*g/1000; %[km/s] exhaust velocity of stage 2
c3 = Is3*g/1000; %[km/s] exhaust velocity of stage 3

c = [c1 c2 c3]'; %effective exhaust velocities vector

%% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DV - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)- c3*log(lambda*c3-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve to find the lagragian multiplier lambda

lambda0 = 1;
lambda = fsolve(fun, lambda0, data.optionsfsolve); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0

%% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage; it is a vector [mr_stage1 mr_stage2 mr_stage3]

m3 = ((mr(3)-1)/(1-eps3.*mr(3)))*m_pay; %[kg] step mass of third stage
m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m3+m_pay); %[kg] step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m3+m2+m_pay); %[kg] step mass of first stage

mi_tot = m1+m2+m3+m_pay; %[kg] total liftoff mass
m = [m1 m2 m3]; %vector of step masses


%% STRUCTURAL MASSES:

ms3 = m3*eps3; %[kg] structural mass of stage 1
ms2 = m2*eps2; %[kg] structural mass of stage 2
ms1 = m1*eps1; %[kg] structural mass of stage 3

ms = [ms1 ms2 ms3]'; %vector of structural masses
%% PROPELLANT MASSES:

mp3 = m3 - ms3;%[kg] propellant mass of stage 1
mp2 = m2 - ms2;%[kg] propellant mass of stage 2
mp1 = m1 - ms1;%[kg] propellant mass of stage 3

mp = [mp1 mp2 mp3]';% vector of propellant masses


%% OUTPUT QUANTITIES
%in this section the masses and indexes computed are saved as outputs as branches of
%the struct launcher, which is a part of the struct data:
launcher.st1.m0 = mi_tot; 
launcher.st1.mp = mp1;
launcher.st1.ms = ms1;

launcher.st2.m0 = m2+m3+m_pay;
launcher.st2.mp = mp2;
launcher.st2.ms = ms2;

launcher.st3.m0 = m3+m_pay;
launcher.st3.mp = mp3;
launcher.st3.ms = ms3;

%% BURNING TIMES AND THRUSTS COMPUTATION:
%in this section the burning time tb is computed for first and second stage and also
%the initial and final thrust the same two stages is found (because the
%values for the third stage have already been fixed into the script
%"config"):
g = orbit.g;%[m/s^2]
stage = fieldnames(launcher);
for i = 1:2 %it applies only on the first and second stage
    
    T_W = launcher.(stage{i}).tw_ratio;
    m0 = launcher.(stage{i}).m0;%[kg]
    mp = launcher.(stage{i}).mp;%[kg]
    isp = launcher.(stage{i}).isp;%[s]
    
    T0 = m0 * T_W * g;%[N]
    Tfinal = (m0-mp) * T_W * g;%[N]
    toll = 1;%first value on the tolerance to be applied at the first iteration of the following while cicle
    mdot_0=T0/(isp*g);%[kg/s]
    mdot_fin = Tfinal/(isp*g);%[kg/s]

    dtb = 0.0001;%discretization on burning time
    tb = 10:dtb:350; %[s]
    j = 0;
    while abs(toll) > 1e-2 && j < length(tb)
        
        j = j + 1;
        mptot = ((mdot_0+mdot_fin)*tb(j))/2;%[kg]
        toll = mp - mptot ;
        
    end
    tb = tb(j);%[s]
     
    launcher.(stage{i}).T = [T0 Tfinal]; %saving the thrust as a vector of initial and final valued in the struct laucher
    launcher.(stage{i}).tb = tb;%saving the burning times in the struct launcher
    
end
%% DeltaVs FOR EACH STAGE:
%here the velocity budgets required by each stage are computed:
DeltaV1 = c1*log(mr(1)); %[m/s] velocity budget required by stage 1
DeltaV2 = c2*log(mr(2)); %[m/s] velocity budget required by stage 2
DeltaV3 = c3*log(mr(3)); %[m/s] velocity budget required by stage 3

DeltaV_stages = [DeltaV1 DeltaV2 DeltaV3]; %vector of velocity budgets
