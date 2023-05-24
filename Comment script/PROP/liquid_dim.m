function [launcher, prop] = liquid_dim(launcher, prop, orbit)

% liquid_dim will provide the sizing of the monopropellant propulsion system and
% will load it in the launcher and prop structs.

% NOZZLE SIZING

eps = prop.eps;     % Area ratio; we want to obtain a nozzle like this 
k = prop.k;         % Heat capacity ratio
eps_inv_t = 1/eps;
p_ratio = 0.00025;  % guess on the pressure ratio
err = 1000;         % initial value of the error, just to start the while

% TRIAL & ERROR procedure to find the correct pressure ratio that satisfies the area ratio chosen
% Starting from the typical expression of the area ratio as fun. of k, p_ratio, using the first guess
% of p_ratio, eps is computed. Then the error is analyzed: if it's greater of +-1% the p_ratio guess
% is decreased/increased till convergence.

while norm(err) > 0.001
     eps_inv = ((k+1)/2)^(1/(k-1))*((p_ratio)^(1/k))*sqrt((k+1)/(k-1)*(1-(p_ratio)^((k-1)/k)));
     err = (eps_inv - eps_inv_t)/ eps_inv_t;
     if err < 0 
         p_ratio = p_ratio + 0.0000005;
     elseif err > 0
         p_ratio = p_ratio - 0.0000005;
     end
end

Pc = prop.Pc;       % Chamber pressure [Pa]
Pe = p_ratio*Pc;    % Exit pressure [Pa]
eps = 1/eps_inv;
    prop.eps = eps;
T0 = prop.T0;       % [K] % Chamber temperature
Ru = prop.Ru;       % Universal gas constant [J/mol*K]
g0 = orbit.g;       % gravity acceleration [m/s^2]
alpha = prop.alpha; % [deg]
beta = prop.beta;   % [deg]
M_c = prop.M_c;     % Mach number in combustion chamber. Assumed value
MM = prop.MM;       % Molar Mass [g/mol]
rho = prop.rho;     % Density [kg/m^3]
T = launcher.T(1);  % Thrust [N]

mprop = launcher.mp;% propellant mass [kg] % safety margin is considered

Pamb = prop.Pamb;   % Ambient pressure. Assumed to 0, third stage is operating in vacuum.

% LRM DESIGN
% thrust coefficient computation
    cf = sqrt(((2*k^2)/(k-1))*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k))) + ((Pe-Pamb)/Pc)*eps;
        stage_3.cf = cf;
    A_t = T/(Pc*cf);          % throat area [m^2]
        stage_3.A_t = A_t;
    r_t = sqrt(A_t/pi);       %throat radius [m]
        stage_3.r_t = r_t;
    A_e = eps*A_t;            % exit area [m^2]
        stage_3.A_e = A_e;
    r_e = sqrt(A_e/pi);       % exit radius [m]
        launcher.r_e = r_e;
%     if 2*r_e > stage_2.dch
%         disp ('ERROR: nozzle cannot fit')
%     end
    
    A_c = A_t/M_c * ((2/(k+1))*(1+(k-1)*0.5*M_c^2))^((k+1)/(2*(k-1))); % Chamber cross section formula [m^2]
        stage_3.A_c = A_c;
    r_c = sqrt(A_c/pi);       % Chamber radius [m]
        stage_3.r_c = r_c;

    v_e = sqrt(2*k/(k-1)*Ru*T0/MM*(1-(Pe/Pc)^((k-1)/k))); % exit velocity % [m/s]
    u_lim = sqrt(2*(k/(k-1))*Ru*T0/MM);                   % limit velocity [m/s]
    eff_noz = 1 - (Pe/Pc)^((k-1)/k);                      % nozzle efficienty
    
    if v_e > u_lim
        disp(' ERROR: exit velocity major than limit velocity')
    elseif v_e == u_lim*eff_noz
        disp('exit velocity ok')
    end
        stage_3.v_e = v_e;
    mdot_p = T/v_e;                % Propellant mass flow rate [kg/s]
        stage_3.mdot_p = mdot_p;
    tburn = mprop/mdot_p;          % Burning time [s]
        launcher.tb = tburn;
    Isp = T/(mdot_p*g0); % [s]     % Specific impulse [s]
    
    toll = 5; % tollerance to check the value of the Isp;
    if norm(Isp - v_e/g0) > toll
        disp ('ERROR: check value of Isp')
    end
    
        stage_3.Isp = Isp;
    Iv = Isp*rho;                  % Volumetric specific impulse 
        launcher.Iv = Iv;
    Itot = T*tburn;                % Total specific impulse
        launcher.Itot = Itot;
    
     % NOZZLE DESIGN: Conical nozzle
    
    lambda =0.5*(1 + cos(alpha));       % Divergence
    
    Ldiv = 0.5*(r_e - r_t)/tan(alpha);  % Divergent Length [m]
        stage_3.Ldiv = Ldiv;
    Lconv = 0.5*(r_c - r_t)/tan(beta);  % Convergent Length [m]
        stage_3.Lconv = Lconv;
    Lnoz = Lconv + Ldiv;                % Nozzle Length [m]
        stage_3.Lnoz = Lnoz;
        
    % Chamber design: we already compute A_c, using the characterstic
    % length. Catalytic bed of platinum-copper bond
    
    L_star = 0.762;           % Characteristic length [m] 
    Vch = L_star*A_t;         % Combustion Chamber volume [m^3]
    Lch = Vch/A_c;            % Combustion chamber length [m]
        stage_3.Lch = Lch;
    
    % Tank size & feeding system design
    
    V_p = mprop/rho;          % Propellant Volume [m^3]
    
    % An increment for temperature variation losses is considered
    
    Vp = V_p*1.03; % 3% increment
     stage_3.Vp = Vp;
    
%     if Vp < 10
%         disp('expectation of  pressure-fed sys compatibility satisfied')
%     end
%     
    % A pressure-fed sys in chosen.
    % In order to compute P tank, losses must be accounted for
    
    injP_loss = Pc*0.1;       % 10% losses for injection [Pa]
    feedP_loss = 50662.5;     % Feeding pipes losses [Pa]
    u_p = 10;                 % Average propellant velocity inside pipes [m/s] % Assumed value    
    dynP_loss = 0.5*rho*u_p^2;% Dynamic losses [Pa]
    
    P_tank = Pc + injP_loss + feedP_loss + dynP_loss; % Tank pressure [Pa]
        stage_3.P_tank = P_tank;
    
    P_burst = 2*P_tank;       % Burst Pressure [Pa] % Safety factor of 2
    
      
    % Pressurising gas (PG) tank design: He
    MM_pg = 4;                                     % PG molar mass [g/mol] 
    k_pg = 1.66;                                   % PG heat capacity ratio
    Ppg_in = 250e+6;                               % PG initial pressure [Pa] % Assumed value
        stage_3.Ppg_in = Ppg_in;
    Ppg_f = P_tank;                                % PG final pressure [Pa] % Pressure is assumed equal to P_tank at end of burning
    Tpg_in = 279.15;                               % PG initial temperature [K]  % guess
    Tpg_f = Tpg_in*(Ppg_f/Ppg_in)^((k_pg-1)/k_pg); % PG final temperature [K] % isoentropic relation 
    Tfreez = 366; % freezing temperature for ADN
    Tboil = 400; % at this temperature ADN start decomposition
%     
%     if Tpg_f < Tfreez
%         disp ('ERROR: gas temperature under freezing point')
%     elseif Tpg_f > Tboil
%         disp ('ERROR: gas temperature may start ADN decomposition')
%     end
%     
    Vpg_f = Vp/(1-(Ppg_f*Tpg_in)/(Ppg_in*Tpg_f)); % PG final volume [m^3] 
    Vpg_in = Vpg_f - Vp;                          % PG initial volume [m^3] % It will be the volume of its tank
        stage_3.Vpg_in = Vpg_in;
    Mpg = Ppg_f*Vpg_f/(Ru*Tpg_f/MM_pg);           % PG mass [kg]
        stage_3.Mpg = Mpg;
    
     
    % for structural considerations: it will now be report the limit value
    % for a diameter, considering all the component in pile (obv not
    % acceptable). this because the third stage must consider also the
    % payload dimensions
    
     Lstage_lim = launcher.C;                  
        launcher.L = Lstage_lim + launcher.fairing;  % modify with toroidal tank design [m]
    
end
    
    
