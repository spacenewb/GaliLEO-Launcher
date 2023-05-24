%%%% STRUCTURAL INDEXES FITTING %%% 
clc
clear 
close all


% stage 1 vega  % vega data sheet 
mp = 87.710;
m = 93.243;
eps(1) = (m-mp)/m;
mpv(1) = mp;

% stage 2 vega
mp = 23.814;
m = 26.300;
eps(2) = (m-mp)/m;
mpv(2) = mp;

% stage 3 vega 
m = 12.000;
mp = 10.567;
eps(3) = (m-mp)/m;
mpv(3) = mp;

%stage1 _pslv http://www.spacelaunchreport.com/pslv.html
m = 168;
mp = 138;
eps(4) = (m-mp)/m;
mpv(4) = mp;


%stage1 _pslv_booster1
m = 10.93;
mp = 8.92;
eps(5) = (m-mp)/m;
mpv(5) = mp;

%stage1 _pslv_booster2
m = 14.7;
mp = 12;
eps(6) = (m-mp)/m;
mpv(6) = mp;

%stage3 _pslv
m = 8.3;
mp = 7.6;
eps(7) = (m-mp)/m;
mpv(7) = mp;

%stage1 _vsb30 https://en.wikipedia.org/wiki/VSB-30
m = 1.200;
mp = 1.200 - 0.341;
eps(8) = (m-mp)/m;
mpv(8) = mp;

%stage2 _vsb30
m = 0.900;
mp = 0.900 - 0.284;
eps(9) = (m-mp)/m;
mpv(9) = mp;


%stage1_shavit9 http://www.spacelaunchreport.com/shavit.html
m = 10.215;
mp = 9.10;
eps(10) = (m-mp)/m;
mpv(10) = mp;

%stage2_shavit9 
m = 10.388;
mp = 9.10;
eps(11) = (m-mp)/m;
mpv(11) = mp;

%stage1_shavit13 
m = 13.99;
mp = 12.75;
eps(12) = (m-mp)/m;
mpv(12) = mp;

%stage2_shavit13 
m = 14.126;
mp = 12.75;
eps(13) = (m-mp)/m;
mpv(13) = mp;

%stage3_shavit 
m = 2.574;
mp = 1.89;
eps(14) = (m-mp)/m;
mpv(14) = mp;

%stage1_mu  https://www.spacelaunchreport.com/m5.html
m = 83.56;
mp = 71.49;
eps(15) = (m-mp)/m;
mpv(15) = mp;

%stage2_mu 
m = 34.47;
mp = 31.06;
eps(16) = (m-mp)/m;
mpv(16) = mp;

%stage3_mu 
m = 11;
mp = 10;
eps(17) = (m-mp)/m;
mpv(17) = mp;

%stage4_mu 
m = 1.43;
mp = 1.31;
eps(18) = (m-mp)/m;
mpv(18) = mp;


%% fitting 
eps = eps';%;[eps(1:4) eps(5:end)]' ;
mpv = mpv';% [mpv(1:4) mpv(5:end)]';
fun = fit(mpv,eps,'poly2');
figure
plot(fun,mpv,eps)
xlabel('m_p [t]')
ylabel('\epsilon')
% p = polyfit(mpv,eps,4);
% eps_v = polyval(p,0:0.1:20);
% 
% figure
% plot(mpv,eps,'o')
% hold on 
% plot(0:0.1:20,eps_v,'-')

