
%% vega

% stage 1 

mp1 = 87.710;
m1 = 93.243;
L = 11.2;
d = 3;
tb1 = 110;
Isp1 = 280;

% stage 2 

mp2 = 23.814;
m2 = 26.300;
L = 8.39;
d = 1.90;
tb2 = 77.1; 
Isp2 = 287.5;

% stage 3 

m3 = 12.000;
mp3 = 10.567;
L = 4.12;
d = 1.90;
tb3 = 119.6; 
Isp3 = 295.9;

% avum 

ms4 = 0.688;
mp4 = 0.381;
m4 = ms4 + mp4;
tb4 = 612.15; 
Isp4 = 314.6;

% calcoli 

eps1 = (m1-mp1)/m1
k1 = (m1-mp1)/mp1;
eps1_t = k1/(1+k1)
eps2 = (m2-mp2)/m2;
eps3 = (m3-mp3)/m3;
eps4 = (m4-mp4)/m4;

mtot = m1+m2+m3+m4+0.540+0.077+1

tw = 3015000/(mtot * 1000 * 9.81)

tw_fin = 3015000/((mtot-mp1) * 1000 * 9.81)


%% Electron 

% stage 1
mp1 = 9.25;
m1 = 10.2;
isp1 = 311;

% stage 2
mp2 = 2.05;
m2 = 2.30;
isp2 = 343;

% satge 3 ???





