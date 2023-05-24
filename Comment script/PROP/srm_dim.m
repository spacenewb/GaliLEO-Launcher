function [launcher,prop] = srm_dim (launcher,prop,orbit,data,index)

% srm_dim provides the sizing of the solid stages and loads their parameters
% into launcher & prop structs. Receives as input also an input (1, 2), which
% will differentiate the nozzle analysis between first and second stage.

T0 = prop.T0;       % Flame temperature [K]
Ru = prop.Ru;       % Universal gas constant [J/mol*K]
g0 = orbit.g;       % Gravity acceleration [m/s^2]
alpha = prop.alpha; % Divergent semiangle [deg]
beta = prop.beta;   % Convergent semiangle [deg]
Pc = prop.Pc;       % Chamber pressure [Pa]
M_c =prop.M_c;      % Combustion chamber Mach number
Tin = launcher.T(1);% Initial THrust value [N]
k = prop.k;         % Heat capacity ratio
MM = prop.MM;       % Molar mass [g/mol]
rho = prop.rho;     % Density [kg/m^3]

% Differentiation between first and second stage

if index == 2               % Second stage: starts solve procedure to find the correct p_ratio that satisfies the area ratio set
    eps = prop.eps;         % Area ratio % set by us
    eps_inv = 1/eps;
    p_ratio0 = 0.00025;     % pressure ratio guess
    fun = @(p_ratio) -eps_inv+((k+1)/2)^(1/(k-1))*((p_ratio)^(1/k))*sqrt((k+1)/(k-1)*(1-(p_ratio)^((k-1)/k)));
    p_ratio = fsolve(fun, p_ratio0, data.optionsfsolve);
    Pe = p_ratio*Pc;        % exit pressure [Pa]
    Pamb = prop.Pamb;       % ambient pressure [Pa]
else                        % First stage procedure (optimal at sea level)
    Pe = prop.Pe;           
    Pamb = Pe;
    eps = 0;                % initialization of area ratio
end


% SRM DESIGN
% thrust coefficient computation
    cf = sqrt(((2*k^2)/(k-1))*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k))) + ((Pe-Pamb)/Pc)*eps;
        stage_1.cf = cf;
    A_t = Tin/(Pc*cf);          % throat area [m^2]
        stage_1.A_t = A_t;
    r_t = sqrt(A_t/pi);         % throat radius [m^2]
        stage_1.r_t = r_t;
    eps_inv = ((k+1)/2)^(1/(k-1))*((Pe/Pc)^(1/k))*sqrt((k+1)/(k-1)*(1-(Pe/Pc)^((k-1)/k))); % A_t / A_e
    eps = 1/eps_inv;
        prop.eps = eps;
    A_e = eps*A_t;              % exit area [m^2]
        stage_1.A_e = A_e;
    r_e = sqrt(A_e/pi);         % exit radius [m]
        launcher.r_e = r_e;
    A_c = A_t/M_c * ((2/(k+1))*(1+(k-1)*0.5*M_c^2))^((k+1)/(2*(k-1))); % Convergent Area [m^2] % this formula is used only for liquid prop. it was used anyway to have a value of the nozzle convergent
        stage_1.A_c = A_c;
    r_c = sqrt(A_c/pi);         % Convergent radius [m]
        stage_1.r_c = r_c;    

    v_e = sqrt(2*k/(k-1)*Ru*T0/MM*(1-(Pe/Pc)^((k-1)/k))); % exit velocity [m/s]
    u_lim = sqrt(2*(k/(k-1))*Ru*T0/MM);                   % limit velocity [m/s]
    eff_noz = 1 - (Pe/Pc)^((k-1)/k);                      % nozzle efficiency
    
    if v_e > u_lim
        disp(' ERROR: exit velocity major than limit velocity')
    elseif v_e == u_lim*eff_noz
        disp('exit velocity ok')
    end
    
    toll = 5; % tollerance to check the value of the Isp;
    if norm(launcher.isp - v_e/g0) > toll
        fprintf ('ERROR: check value of Isp %d \n',index)
    end
    
    Iv = launcher.isp*rho;      % Volumetric specific mpulse
        launcher.Iv = Iv;
    
    
    % NOZZLE DESIGN (first stage): Bell-shaped nozzle to reduce divergence losses
    % for eps around 9, choosing max performance from RAO charts:
    
    th_i = 21.9;            % Initial inclination  [deg]
    th_i = deg2rad(th_i);
    th_e = 6.6;             % Final inclination [deg]
    th_e = deg2rad(th_e);
    lambda =0.5*(1 + cos(0.5*(alpha + th_e))); % Divergence
    
    Ldiv = 0.5*(r_e - r_t)/tan(alpha);  % Divergent length [m]
        stage_1.Ldiv = Ldiv;
    Lconv = 0.5*(r_c - r_t)/tan(beta);  % Convergent length [m]
        stage_1.Lconv = Lconv;
    Lnoz = Lconv + Ldiv;                % Nozzle length [m]
        launcher.Lnoz = Lnoz;
    
    
    
    % GRAIN CONFIGURATION
    
    
    V_p = launcher.mp/rho;         % Propellant Volume
        stage_1.V_p = V_p;
    a = 1.21;                      % Vieille law constant
    n = 0.41;                      % Vieille law pressure exponent
    r_b = a*(Pc*10^-5)^n;          % BUrning rate [mm/s] Pressure in bar
        stage_1.r_b = r_b;
        
    % we consider now a Double Anchor configuration for a regressive behaviour: from tables we can find
    % the volumetric fraction
    
    V_f = 0.725;            % Volumetric fraction
    Vch = V_p / V_f;        % Volume of the combustion chamber [m^3]
        stage_1.Vch = Vch;
   
    % A guess on the diameter has been made and considering a range of acceptable L/D
    % the length of the combustion chamber is computed
    
    dch = launcher.C;           % Chamber diameter [m] assumption
    Lch = Vch/((dch/2)^2 * pi); % Combustion chamber length [m]
        launcher.Lcc = Lch;
        
end

    
    
    
    
    
    
    
