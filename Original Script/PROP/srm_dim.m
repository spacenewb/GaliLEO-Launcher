function [launcher,prop] = srm_dim (launcher,prop,orbit,data,index)


T0 = prop.T0; % [K]
Ru = prop.Ru;
g0 = orbit.g;
alpha = prop.alpha; % [deg]
beta = prop.beta; % [deg]
M_c = prop.M_c;
Pc = prop.Pc;
Tin = launcher.T(1);
k = prop.k;
MM = prop.MM;
rho = prop.rho;


if index == 2
    eps = prop.eps;
    eps_inv = 1/eps;
    p_ratio0 = 0.00025;
    fun = @(p_ratio) -eps_inv+((k+1)/2)^(1/(k-1))*((p_ratio)^(1/k))*sqrt((k+1)/(k-1)*(1-(p_ratio)^((k-1)/k)));
    p_ratio = fsolve(fun, p_ratio0, data.optionsfsolve);
    Pe = p_ratio*Pc;
    Pamb = prop.Pamb; %Pa
else 
    Pe = prop.Pe;
    Pamb = Pe;
    eps = 0;
end


% SRM DESIGN
% thrust coefficient computation
    cf = sqrt(((2*k^2)/(k-1))*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k))) + ((Pe-Pamb)/Pc)*eps;
        stage_1.cf = cf;
    A_t = Tin/(Pc*cf); % throat area [m^2]
        stage_1.A_t = A_t;
    r_t = sqrt(A_t/pi);
        stage_1.r_t = r_t;
    eps_inv = ((k+1)/2)^(1/(k-1))*((Pe/Pc)^(1/k))*sqrt((k+1)/(k-1)*(1-(Pe/Pc)^((k-1)/k))); % A_t / A_e
    eps = 1/eps_inv;
        prop.eps = eps;
    A_e = eps*A_t; % exit area [m^2]
        stage_1.A_e = A_e;
    r_e = sqrt(A_e/pi);
        launcher.r_e = r_e;
    A_c = A_t/M_c * ((2/(k+1))*(1+(k-1)*0.5*M_c^2))^((k+1)/(2*(k-1)));
        stage_1.A_c = A_c;
    r_c = sqrt(A_c/pi);
        stage_1.r_c = r_c;    

    v_e = sqrt(2*k/(k-1)*Ru*T0/MM*(1-(Pe/Pc)^((k-1)/k))); % [m/s]
    u_lim = sqrt(2*(k/(k-1))*Ru*T0/MM);
    eff_noz = 1 - (Pe/Pc)^((k-1)/k);
    
    if v_e > u_lim
        disp(' ERROR: exit velocity major than limit velocity')
    elseif v_e == u_lim*eff_noz
        disp('exit velocity ok')
    end
    
    toll = 5; % tollerance to check the value of the Isp;
    if norm(launcher.isp - v_e/g0) > toll
        fprintf ('ERROR: check value of Isp %d \n',index)
    end
    
    Iv = launcher.isp*rho;
        launcher.Iv = Iv;
    
    
    % NOZZLE DESIGN: Bell-shaped nozzle to reduce divergence losses
    % for eps around 9, choosing max performance:
    
    th_i = 21.9; % [deg]
    th_i = deg2rad(th_i);
    th_e = 6.6; % [deg]
    th_e = deg2rad(th_e);
    lambda =0.5*(1 + cos(0.5*(alpha + th_e)));
    
    Ldiv = 0.5*(r_e - r_t)/tan(alpha);
        stage_1.Ldiv = Ldiv;
    Lconv = 0.5*(r_c - r_t)/tan(beta);
        stage_1.Lconv = Lconv;
    Lnoz = Lconv + Ldiv;
        launcher.Lnoz = Lnoz;
    
    
    
    % GRAIN CONFIGURATION
    
    
    V_p = launcher.mp/rho;
        stage_1.V_p = V_p;
    a = 1.21;
    n = 0.41;
    r_b = a*(Pc*10^-5)^n;% [mm/s] Pressure in bar
        stage_1.r_b = r_b;
        
    % we consider now a Double Anchor configuration for a regressive behaviour: from tables we can find
    % the volumetric fraction
    
    V_f = 0.725;
    Vch = V_p / V_f; % volume of the combustion chamber
        stage_1.Vch = Vch;
   
    % a guess on the length of the first stage has to be made in order to
    % determine the internal diameter of the case. a comparison with Vega
    % is considered in the choice of L/d
    
    dch = launcher.C; % [m]
    Lch = Vch/((dch/2)^2 * pi); % [m]
        launcher.Lcc = Lch;
        
end

    
    
    
    
    
    
    
