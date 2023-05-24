clc
clear all
close all

run config.m

% Recovery --> 1 main parachute + 1 or more drogue pilot chutes
%% Data

g0 = 9.80665; % m/s^2
R_Earth = 6371; % km
m_payload = 250 ; % Mass of 1st stage
w_payload = m_payload*g0 ; % Weight of first stage

% Touch down data (constrains)

v_td = 5; % [m/s] Touch down velocity
z_td = 0; % [m] Altitude wrt sea level

% Opening data

v_open = 2115; % Speed at opening of 1st stage
z_open = 55526; % Opening altitude 
g_open = g0*(R_Earth/(R_Earth+z_open*1e-3))^2; % Opening gravity
a = 6;      % Shock resistance 
[~,~,~,rho_open] = atmoscoesa(z_open,'none');
[~,~,~,rho0] = atmoscoesa(0);
v_eas_open = v_open * sqrt(rho_open/rho0);
v_eas_Kn = v_eas_open * 1.944;

%% Main parachute

c_d_main = 0.8; % Conical shape parachute
c_x_main = 1.8; % Opening force coefficient (for conical shape)

[~,~,~,rho_td] = atmoscoesa(z_td,'None'); % Density at sea level (touch down)

q_td = .5*rho_td*v_td^2; % Dynamic pressure
D = w_payload; % Drag in equilibrium condition at opening = weight of payload
S_main = D/(q_td*c_d_main); % Parachute surface
d_main = sqrt((4*S_main)/pi); % Parachute diameter

W_main = 10; % Took out from table diameter-weight on Parachute slide n°9. 
             % If W > 2kg we need one ore more pilot chutes
             
%% Drogue of pilot parachute (if W>2 we need another one)

U = 0.03; % Coefficient for pilot to main parachute area ratio for our speed
S_drogue = 4*U*S_main;
d_drogue = sqrt((4*S_drogue)/pi);

% From these values we can find the type of drogue pilot desing (from the
% diameter and/or surface) and define:
c_d_drogue = 0.9;
c_x_drogue = 1.05;

W_drogue = 1.8; % From tables


% % WE HAVE TO CHECK IF W_DROGUE < 2, IF NOT WE NEED ANOTHER PILOT CHUTE SO WE HAVE
% % TO COMPUTE THE SAME PARAMETERS TROUGH THE SAME PASSAGES
% % IF WE NEED IT
% 
% % Drogue pilot chute design 
% % 
% % c_d_drogue_p = ??;
% % c_x_drogue_p = ??;
% % 
% % S_drogue_p = U*X*S_drogue;  %pilot to main paracute area ratio for high speed = 0.01
% % d_drogue_p = sqrt((4*S_drogue_p)/pi);
% % 
% % W_drogue1 = ??;
% 
% 
% %% Main opening altitude estimation + Pflanz method for dynamic forces on parachute
% 
% S_0_drogue = //;         % D_c/D_0 = ?? from table 
% S_r = S_main;            % For unreefed parachute;
% n = //;                  % For unreefed conical ribbon
% 
% W_t = w_payload;
% z = linspace(0, z_open, 1000);
% 
% for i=1:length(z)
%     [~,~,~,rho_main(i)] = atmoscoesa(z(i),'None');
%     g_m(i) = g0*(6371/(6371+z(i)*1e-3))^2;
%     v_e(i) = sqrt((2*W_t)/(c_d_drogue*S_0_drogue*rho_main(i)));   % Rate of descent
%     t(i) = ((n*d_main)/v_e(i))*sqrt(S_r/S_main);                  % Filling time
%     A(i) = (2*W_t)/(S_r*rho_main(i)*g_m(i)*t(i)*v_e(i));          % Ballistic factor   
% end
% 
% X = ??; % Scaling factor check on tables    [???????????]
% F = 0.5*rho_main(68)*v_e(68)^2*S_main*c_x_main*X; % Pflanz method
% a = F/w_payload; % Shock
% 
% %% Drogue opening altitude estimation + Pflanz method for dynamic forces on pilot chute
% 
% S_0_drogue_p = S_drogue/1;   % D_c/D_0 = ?? from table    ----- S_drogue_p is referred to the pilot chute
% S_r_d = S_drogue;            % For unreefed parachute;
% n_d = //;                    % For unreefed conical ribbon
% 
% W_t_d = w_payload + W_main;
% 
% for i=1:length(z)
%     [~,~,~,rho_m_d(i)] = atmoscoesa(z(i),'None');
%     g_m_d(i) = g0*(6371/(6371+z(i)*1e-3))^2;
%     v_e_d(i) = sqrt((2*W_t_d)/(c_d_drogue*S_0_drogue_p*rho_m_d(i))); 
%     t_d(i) = ((n_d*d_drogue)/v_e_d(i))*sqrt(S_r_d/S_drogue);         
%     A_d(i) = (2*W_t_d)/(S_r_d*rho_m_d(i)*g_m_d(i)*t_d(i)*v_e_d(i));  
% end
% 
% X_d = ??; % Scaling factor for chutes  --> We have to check the opening altitude
% F_d = 0.5*rho_m_d(278)*v_e_d(278)^2*S_drogue*c_x_drogue*X_d;
% a_d = F_d/w_payload
% 
% S_r_p = S_drogue_p;
% t_d_p = ((n_d*d_drogue_p)/v_open)*sqrt(S_r_p/S_drogue_p);
% 
% % From these we have to took out the altitude of pilot and main chute,
% % velocities, weights, t, A, X, F[N], a (shock)
% 
