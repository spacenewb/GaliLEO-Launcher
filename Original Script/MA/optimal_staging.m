function [launcher,DeltaV_stages] = optimal_staging(data)

orbit = data.orbit;
launcher = data.launcher;

%% DATA:

rp = orbit.rp_e + orbit.Re;
ra = orbit.h + orbit.Re;
e = (ra-rp)/(ra+rp);
a = (ra+rp)/2;
p = a*(1-e^2);
vp = sqrt(orbit.mu/p)*(1+e);
va = sqrt(orbit.mu/p)*(1-e);

v_circ = sqrt( orbit.mu/(orbit.h+orbit.Re) ) ;
dv_orbit = v_circ - va;


dv_loss = orbit.losses/1000;
DV = (vp + dv_orbit + dv_loss);
g = orbit.g; %[m/s] GRAVITY ACCELERATION ON EARTH

m_pay = launcher.st3.pay; %[kg] PAYLOAD MASS

%structural coefficients:
eps1 = launcher.st1.eps;
eps2 = launcher.st2.eps;
eps3 = launcher.st3.eps;

eps = [eps1 eps2 eps3]'; %structural coefficients vector

%Specific impulses:
Is1 = launcher.st1.isp; %[s]
Is2 = launcher.st2.isp; %[s]
Is3 = launcher.st3.isp; %[s]

%Exhaust velocities:
c1 = Is1*g/1000; %[km/s]
c2 = Is2*g/1000; %[km/s]
c3 = Is3*g/1000; %[km/s]

c = [c1 c2 c3]'; %effective exhaust velocities vector

%% LAGRAGIAN MULTIPLIERS COMPUTATION:

fun = @(lambda) DV - c1*log(lambda*c1-1)- c2*log(lambda*c2-1)- c3*log(lambda*c3-1)+log(lambda)*sum(c)+sum(c.*log(c.*eps)); %function which we have to solve
%fun = @(lambda) DV - sum(c.*log( (lambda*c - 1)./(lambda.*c.*eps) ));

lambda0 = 1;
lambda = fsolve(fun, lambda0, data.optionsfsolve); %solves fun with the initial condition lambda0 ATTENTION: very sensitive to lambda0

%% STEP MASSES & INITIAL TOTAL MASS COMPUTATION:

mr = (lambda.*c-1)./(lambda.*c.*eps); %mass ratios for each stage

m3 = ((mr(3)-1)/(1-eps3.*mr(3)))*m_pay; %step mass of third stage
m2 =((mr(2)-1)/(1-eps2.*mr(2)))*(m3+m_pay); %step mass of second stage
m1 =((mr(1)-1)/(1-eps1.*mr(1)))*(m3+m2+m_pay); %step mass of first stage

mi_tot = m1+m2+m3+m_pay;
m = [m1 m2 m3];


%% STRUCTURAL MASSES:

ms3 = m3*eps3;
ms2 = m2*eps2;
ms1 = m1*eps1;

ms = [ms1 ms2 ms3]';
%% PROPELLANT MASSES:

mp3 = m3 - ms3;
mp2 = m2 - ms2;
mp1 = m1 - ms1;

mp = [mp1 mp2 mp3]';


%% OUTPUT QUANTITIES
launcher.st1.m0 = mi_tot;
launcher.st1.mp = mp1;
launcher.st1.ms = ms1;

launcher.st2.m0 = m2+m3+m_pay;
launcher.st2.mp = mp2;
launcher.st2.ms = ms2;

g = orbit.g;
stage = fieldnames(launcher);
for i = 1:2
    
    T_W = launcher.(stage{i}).tw_ratio;
    m0 = launcher.(stage{i}).m0;
    mp = launcher.(stage{i}).mp;
    isp = launcher.(stage{i}).isp;
    
    T0 = m0 * T_W * g;
    Tfinal = (m0-mp) * T_W * g;
    toll = 1;
    mdot_0=T0/(isp*g);
    mdot_fin = Tfinal/(isp*g);

    dtb = 0.0001;
    tb = 10:dtb:350;
    j = 0;
    while abs(toll) > 1e-2 && j < length(tb)
        
        j = j + 1;
        mptot = ((mdot_0+mdot_fin)*tb(j))/2;
        toll = mp - mptot ;
        
    end
    tb = tb(j);
     
    launcher.(stage{i}).T = [T0 Tfinal];
    launcher.(stage{i}).tb = tb;
    
end

launcher.st3.m0 = m3+m_pay;
launcher.st3.mp = mp3;
launcher.st3.ms = ms3;

%% DeltaVs FOR EACH STAGE:

DeltaV1 = c1*log(mr(1));
DeltaV2 = c2*log(mr(2));
DeltaV3 = c3*log(mr(3));

DeltaV_stages = [DeltaV1 DeltaV2 DeltaV3];
