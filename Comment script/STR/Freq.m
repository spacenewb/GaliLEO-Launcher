function Freq(launcher,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates the evolution of the first modal frequency throughout the mission 
% duration. It also calls the function responsible to evaluate the modal frequency and shape
% of the launcher using Free-Free beam model, at the initial launcher configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Value and Constant Declarations

M0 = launcher.st1.m0 + launcher.st2.m0 + launcher.st3.m0 + launcher.st3.pay;    % Total mass of launcher

m_dot = [(launcher.st1.mp/launcher.st1.tb) (launcher.st2.mp/launcher.st2.tb) (launcher.st3.mp/launcher.st3.tb)];    % Mass flow rate vector

tbb = [launcher.st1.tb launcher.st2.tb launcher.st3.tb];    % Time of burn vector

L = [14.52 1.95 1.83];  % Stage length vector

Diam = 0.83;    % Diameter of launcher

thick = 0.005;  % Thickness of external shell

E = data.str.E;     % Youngs Modulus of the material of external shell


%% Initialising variable values based on the stage that is operational at a given time of burn

mass=zeros(floor(sum(tbb)),1);
time=zeros(floor(sum(tbb)),1);
Length=zeros(floor(sum(tbb)),1);

for mission_time=0:sum(tbb)
    
    if mission_time < tbb(1)
        
        mass(mission_time+1) = M0 - (m_dot(1)*mission_time);
        
        time(mission_time+1) = mission_time;
        
        Length(mission_time+1) = sum(L);
        
    elseif (mission_time >= tbb(1)) && (mission_time < sum(tbb(1:2)))
        
        mass(mission_time+1) = M0 - launcher.st1.m0 - (m_dot(2)*(mission_time - tbb(1))); 
        
        time(mission_time+1) = mission_time;
        
        Length(mission_time+1) = sum(L(2:end));
        
    else
        
        mass(mission_time+1) = M0 - launcher.st1.m0 - launcher.st2.m0 - (m_dot(3)*(mission_time - tbb(1) - tbb(2)));
        
        time(mission_time+1) = mission_time;
        
        Length(mission_time+1) = sum(L(end));
        
    end
    
end

%% evaluations

LoD_time = Length./Diam;    % Current Fineness ratio at a given time during the mission

Et_M = (E.*thick)./mass;    

f1 = (9.87./(2.*pi)).*sqrt((pi.*Et_M)./(8.*(LoD_time.^3))); % Current First modal frequency at a given time during the mission

Rho_eq = 5*(4*M0)/(pi*sum(L)*(((Diam+(2*thick))^2) - (Diam^2)));    % Current equivalent density (using weighted average) at a given time during the mission

how_many_modes = 2; % Number of modes to calculate

[xl,Xnx,~] = FFbeam(Diam,sum(L),thick,E,Rho_eq,how_many_modes); % Calling the function to evaluate the modal frequencies and shapes



%% plots

mission_burn_time = time;

figure
    plot(mission_burn_time,f1,'r','Linewidth',1)
    title('Mode 1 Frequency vs Burn Time')
    xlabel('Burn Time [s]');
    ylabel('First mode frequency [Hz]');
    grid on
    

figure
    img = imread('Rocket_sketch.png');
    image('CData',img,'XData',[-(Diam/2/sum(L)) (Diam/2/sum(L))],'YData',[sum(L)/sum(L) 0])
    hold on
    plot(Xnx(1,:),xl, 'b-',Xnx(2,:),xl, 'r-',[0 0],[0 1.1], '--k');
    title('Bending Mode Shapes')
    subtitle('(Free-Free Model)')
    xlabel('Mode Shape [ ]');
    ylabel('Normalised Length [ ]');
    axis equal
    ylim([0 1.1])
    %xlim([-0.1 0.1])
    legend('Mode 1','Mode 2','Location','northeast');
    grid on
    hold off


end
