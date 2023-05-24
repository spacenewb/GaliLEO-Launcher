function [gamma_obj] = gamma_fit(dt)
state = load('orbit_path.mat');
Y = state.Y;
T = state.T;
gamma = Y(:,4);

%% fitting

% f = fit(T,gamma,ft,'problem',4);
% gamma_obj = feval(f,T(1):dt:T(end));

% p = polyfit(T,gamma,16);
% gamma_obj = polyval(p,T(1):dt:T(end));

 gamma_obj = csaps(T,gamma,0.1,T(1):dt:T(end));
 

