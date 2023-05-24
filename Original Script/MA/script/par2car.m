%funzione pasaggio da parametriche a cartesiane 

function [rr,vv]=par2car(a,e,i,OM,om,th,mu)


ii=[1,0,0];
jj=[0,1,0];
kk=[0,0,1];

% passo 1 calcolo r in perifocale 

r=a*(1-e^2)/(1+e*cos(th));

rrp=[r*cos(th),r*sin(th),0];


% passo 2 calcolo v in perifocale 

p=a*(1-e^2);
vvp=[-sin(th) e+cos(th) 0].*sqrt(mu/p);

% passo 3 rperifocale-geocentrico inerziale 

% definisco matrici di rotazione

ROM=[cos(OM) sin(OM) 0
    -sin(OM) cos(OM) 0
    0 0 1];
Ri=[1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
Rom=[cos(om) sin(om) 0 
    -sin(om) cos(om) 0
    0 0 1];

T=Rom*Ri*ROM;

% ruoto r 

rr=T'*rrp';

% ruoto v

vv=T'*vvp';
end





