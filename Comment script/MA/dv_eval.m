function cost = dv_eval(x,launcher,orbit,site,N,options)

% setting initial condition of integration
x0 = 0;
z0 = site.z0;
V0 = 1;
gamma0 = deg2rad(site.OMEGA);

% set time burning of 3rd stage 
launcher.st3.tb = x(1);
% set control
chi = x(2:end);

% initialize variable for the cycle 
stage = fieldnames(launcher); 
k = 0;
t0_stage = 0;
Y_tot = []; 
T_tot = [];

% cycle among stages 
for i = 1:numel(stage)
    
    m0 = launcher.(stage{i}).m0;                    % set initila mass for actual stage 
    Y0 = [x0 z0 V0 gamma0 m0];                      % initial state vector 
    tf = launcher.(stage{i}).tb + t0_stage;         % update final time for ode
    launcher.(stage{i}).T_time = [t0_stage tf];     % set the burning tme for actual stage 
    
    dt = (tf-t0_stage)/N;

% cycle for N control variable 
for j = 1:N 
   
    t0 = (j-1)*dt + t0_stage; 
    tf = j*dt + t0_stage;
    [T, Y] = ode45(@ascent_opt, [t0 tf], Y0, options, launcher.(stage{i}), site, orbit, chi(j+k));
    Y0 = Y(end,:); % update initial state
    % assemble state 
    Y_tot = [Y_tot; Y];
    T_tot = [T_tot; T];

end

    k = j;                      % counter for control vector
    % update initial condition
    x0 = Y(end,1);
    z0 = Y(end,2);
    V0 = Y(end,3);
    gamma0 = Y(end,4);
    t0_stage = tf; 
 
   
end

   
    CD = 0.49;      % CD avarage (evaluated with base simulator for a similar trajectory)
    % evaluation of the losses 
    [~,~,~,rho] = atmoscoesa(Y_tot(:,2),'none');
    V = Y_tot(:,3);
    Re = orbit.Re * 1e3; % Radius of Earth
    drag = 0.5*rho.*V.^2.*launcher.st1.S.*CD;
    drag_loss = drag./Y_tot(:,5);                                           
    index = not(isnan(drag_loss));
    dv_drag = trapz(T_tot(index),drag_loss(index));                         % drag losses 
    grav_loss = orbit.g .* ( Re./(Re + Y_tot(:,2)) ).^2 .* sin(Y_tot(:,4)); 
    dv_grav = trapz(T_tot,grav_loss);                                       % grav losses
    
    % evaluate cost 
    cost = dv_grav+dv_drag; % [m/s]




