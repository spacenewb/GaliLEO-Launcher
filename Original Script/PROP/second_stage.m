function [stage_2] = second_stage(T, Tfin, Pc, MM, k, rho, mprop, stage_1)

% second stage will operate in vacuum, so nozzle is not working in optimal
% expansion. A trial  & error procedure is adopted to find the pressure
% ratio, starting with a plausible guess of epsilon and a first guess of
% the pressure ratio

eps = 40; % we want to obtain a nozzle like this
eps_inv_t = 1/eps;
p_ratio = 0.00025;
k = 1.1841;
err = 1000;

while norm(err) > 0.001
     eps_inv = ((k+1)/2)^(1/(k-1))*((p_ratio)^(1/k))*sqrt((k+1)/(k-1)*(1-(p_ratio)^((k-1)/k)));
     err = (eps_inv - eps_inv_t)/ eps_inv_t;
     if err < 0 
         p_ratio = p_ratio + 0.0000005;
     elseif err > 0
         p_ratio = p_ratio - 0.0000005;
     end
end

Pe = p_ratio*Pc;
eps = 1/eps_inv;
    stage_2.eps = eps;

T0 = 2643.3; % [K]
Ru = 8314.46;
g0 = 9.81;
alpha = 15; % [deg]
beta = 45; % [deg]
alpha = deg2rad(alpha);
beta = deg2rad(beta);
M_c = 0.2;

m_prop = 1.1*mprop; % propellant mass increment for robust design

% we want the nozzle to be optimal at sea level so Pe = Pa at z = 0
Pamb = 0; %Pa

% SRM DESIGN
% thrust coefficient computation
     stage_2.T = T;
    cf = sqrt(((2*k^2)/(k-1))*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k))) + ((Pe-Pamb)/Pc)*eps;
        stage_2.cf = cf;
    A_t = T/(Pc*cf); % throat area [m^2]
        stage_2.A_t = A_t;
    r_t = sqrt(A_t/pi);
        stage_2.r_t = r_t;
    A_e = eps*A_t; % exit area [m^2]
        stage_2.A_e = A_e;
    r_e = sqrt(A_e/pi);
        stage_2.r_e = r_e;
    if 2*r_e > stage_1.dch
        disp ('ERROR: nozzle cannot fit')
    end
    
    A_c = A_t/M_c * ((2/(k+1))*(1+(k-1)*0.5*M_c^2))^((k+1)/(2*(k-1)));
        stage_2.A_c = A_c;
    r_c = sqrt(A_c/pi);
        stage_2.r_c = r_c;

    v_e = sqrt(2*k/(k-1)*Ru*T0/MM*(1-(Pe/Pc)^((k-1)/k))); % [m/s]
    u_lim = sqrt(2*(k/(k-1))*Ru*T0/MM);
    eff_noz = 1 - (Pe/Pc)^((k-1)/k);
    
    if v_e > u_lim
        disp(' ERROR: exit velocity major than limit velocity')
    elseif v_e == u_lim*eff_noz
        disp('exit velocity ok')
    end
        stage_2.v_e = v_e;
    mdot_p = T/v_e; % [kg/s]
        stage_2.mdot_p = mdot_p;
     mdot_p_f = Tfin/v_e; % [kg/s]
        stage_2.mdot_p_f = mdot_p_f;
    tburn = m_prop/mdot_p;
        stage_2.tburn = tburn;
    Isp = T/(mdot_p*g0); % [s]
    
    toll = 5; % tollerance to check the value of the Isp;
    if norm(Isp - v_e/g0) > toll
        disp ('ERROR: check value of Isp')
    end
    
        stage_2.Isp = Isp;
    Iv = Isp*rho;
        stage_2.Iv = Iv;
    Itot = T*tburn;
        stage_2.Itot = Itot;
    
    
    % NOZZLE DESIGN: Conical nozzle
    
    lambda =0.5*(1 + cos(alpha));
    
    Ldiv = 0.5*(r_e - r_t)/tan(alpha);
        stage_2.Ldiv = Ldiv;
    Lconv = 0.5*(r_c - r_t)/tan(beta);
        stage_2.Lconv = Lconv;
    Lnoz = Lconv + Ldiv;
        stage_2.Lnoz = Lnoz;
    
    
    
    % GRAIN CONFIGURATION
    
    V_p = m_prop/rho;
        stage_2.V_p = V_p;
    a = 1.21;
    n = 0.41;
    r_b = a*(Pc*10^-5)^n;% [mm/s] Pressure in bar
        stage_2.r_b = r_b;
    rb = r_b*10^-3; % coversion from mm/s to m/s
    A_b = mdot_p/(rho*rb); % [m^2]
        stage_2.A_b = A_b;
     A_b_f = mdot_p_f/(rho*rb);
        stage_2.A_b_f = A_b_f;
            
    % we consider now a Double Anchor configuration for a regressive behaviour: from tables we can find
    % the volumetric fraction
    
    V_f = 0.725;
    Vch = V_p / V_f; % volume of the combustion chamber
   
    % a guess on the length of the second stage stage has to be made in order to
    % determine the internal diameter of the case. a comparison with Vega
    % is considered in the choice of L/d
    
    fin_ratio = 5;
        stage_2.fin_ratio = fin_ratio;
    dch = (Vch*4/(fin_ratio*pi))^(1/3); % [m]
        stage_2.dch = dch;
    Lch = fin_ratio*dch; % [m]
        stage_2.Lch = Lch;
    
    Klem = A_b/A_t;

if Klem < 5
    disp('ATTENTION: erosive burning');
end

        
end
