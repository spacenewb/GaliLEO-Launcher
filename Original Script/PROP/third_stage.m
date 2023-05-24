function [stage_3] = third_stage(T, Pc, MM, k, rho, mprop)

eps = 200; % we want to obtain a nozzle like this
eps_inv_t = 1/eps;
p_ratio = 0.00025;
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
    stage_3.eps = eps;
T0 = 2131.11; % [K]
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

% LRM DESIGN
% thrust coefficient computation
     stage_3.T = T;
    cf = sqrt(((2*k^2)/(k-1))*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k))) + ((Pe-Pamb)/Pc)*eps;
        stage_3.cf = cf;
    A_t = T/(Pc*cf); % throat area [m^2]
        stage_3.A_t = A_t;
    r_t = sqrt(A_t/pi);
        stage_3.r_t = r_t;
    A_e = eps*A_t; % exit area [m^2]
        stage_3.A_e = A_e;
    r_e = sqrt(A_e/pi);
        stage_3.r_e = r_e;
%     if 2*r_e > stage_2.dch
%         disp ('ERROR: nozzle cannot fit')
%     end
    
    A_c = A_t/M_c * ((2/(k+1))*(1+(k-1)*0.5*M_c^2))^((k+1)/(2*(k-1)));
        stage_3.A_c = A_c;
    r_c = sqrt(A_c/pi);
        stage_3.r_c = r_c;

    v_e = sqrt(2*k/(k-1)*Ru*T0/MM*(1-(Pe/Pc)^((k-1)/k))); % [m/s]
    u_lim = sqrt(2*(k/(k-1))*Ru*T0/MM);
    eff_noz = 1 - (Pe/Pc)^((k-1)/k);
    
    if v_e > u_lim
        disp(' ERROR: exit velocity major than limit velocity')
    elseif v_e == u_lim*eff_noz
        disp('exit velocity ok')
    end
        stage_3.v_e = v_e;
    mdot_p = T/v_e; % [kg/s]
        stage_3.mdot_p = mdot_p;
    tburn = m_prop/mdot_p;
        stage_3.tburn = tburn;
    Isp = T/(mdot_p*g0); % [s]
    
    toll = 5; % tollerance to check the value of the Isp;
    if norm(Isp - v_e/g0) > toll
        disp ('ERROR: check value of Isp')
    end
    
        stage_3.Isp = Isp;
    Iv = Isp*rho;
        stage_3.Iv = Iv;
    Itot = T*tburn;
        stage_3.Itot = Itot;
    
     % NOZZLE DESIGN: Conical nozzle
    
    lambda =0.5*(1 + cos(alpha));
    
    Ldiv = 0.5*(r_e - r_t)/tan(alpha);
        stage_3.Ldiv = Ldiv;
    Lconv = 0.5*(r_c - r_t)/tan(beta);
        stage_3.Lconv = Lconv;
    Lnoz = Lconv + Ldiv;
        stage_3.Lnoz = Lnoz;
        
    % Chamber design: we already compute A_c, using the characterstic
    % length. Catalytic bed of platinum-copper bond
    L_star = 0.762;
    Vch = L_star*A_t;
    Lch = Vch/A_c;
        stage_3.Lch = Lch;
    
    % Tank size & feeding system design
    
    V_p = mprop/rho;
    
    % An increment for temperature variation losses is considered
    
    Vp = V_p*1.03; % 3% increment
     stage_3.Vp = Vp;
    
    if Vp < 10
        disp('expectation of  pressure-fed sys compatibility satisfied')
    end
    
    % A pressure-fed sys in chosen.
    % In order to compute P tank, losses must be accounted for
    
    injP_loss = Pc*0.1; % 10% losses for injection
    feedP_loss = 50662.5; % [Pa]
    u_p = 10; % [m/s]
    dynP_loss = 0.5*rho*u_p^2; % [Pa]
    
    P_tank = Pc + injP_loss + feedP_loss + dynP_loss; % [Pa]
        stage_3.P_tank = P_tank;
    
    % we choose spherical tanks
    
    r_tank = (3*Vp/(4*pi))^(1/3); 
    A_tank = 4*pi*r_tank^2; 
    
    P_burst = 2*P_tank; % safety factor of 2
    
    % We choose tanks made of carbon fiber since it's the most structurally
    % efficient
    rho_m = 1550;
    Ftu_m = 0.895e+9;
    t_tank = P_burst*r_tank/(2*Ftu_m);
    M_tank = A_tank.*t_tank.*rho_m;
        stage_3.M_tank = M_tank;
    
    % pressurising gas tank design: He
    MM_pg = 4;
    k_pg = 1.66;
    Ppg_in = 250e+6;
        stage_3.Ppg_in = Ppg_in;
    Ppg_f = P_tank;
    Tpg_in = 279.15; % [K] guess
    Tpg_f = Tpg_in*(Ppg_f/Ppg_in)^((k_pg-1)/k_pg);
    Tfreez = 366; % freezing temperature for ADN
    Tboil = 400; % at this temperature ADN start decomposition
    
    if Tpg_f < Tfreez
        disp ('ERROR: gas temperature under freezing point')
    elseif Tpg_f > Tboil
        disp ('ERROR: gas temperature may start ADN decomposition')
    end
    
    Vpg_f = Vp/(1-(Ppg_f*Tpg_in)/(Ppg_in*Tpg_f));
    Vpg_in = Vpg_f - Vp;
        stage_3.Vpg_in = Vpg_in;
    Mpg = Ppg_f*Vpg_f/(Ru*Tpg_f/MM_pg);
        stage_3.Mpg = Mpg;
    
    % Pressurizing tank design (spherical, made of carbon fiber as before)
    rpg_tank = (3*Vpg_in/(4*pi))^(1/3);
    Apg_tank = 4*pi*rpg_tank^2;
    tpg_tank = 2*Ppg_in*rpg_tank/(2*Ftu_m);
    Mpg_tank = Apg_tank*rho_m*tpg_tank;
        stage_3.Mpg_tank = Mpg_tank;
    
    % for structural considerations: it will now be report the limit value
    % for a diameter, considering all the component in pile (obv not
    % acceptable). this because the third stage must consider also the
    % payload dimensions
    diam = [2*r_c; 2*r_tank; 2*rpg_tank];
    dstage_lim = max(diam);
        stage_3.dstage_lim = dstage_lim;
     Lstage_lim = Lch + 2*r_tank + 2*rpg_tank;
        stage_3.Lstage_lim = Lstage_lim;  % modify with toroidal tank design 
    
end
    
    
