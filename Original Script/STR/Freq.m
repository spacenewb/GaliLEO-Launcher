function Freq(launcher,data)

%launcher = data.launcher;
%%

M0 = launcher.st1.m0 + launcher.st2.m0 + launcher.st3.m0 + launcher.st3.pay;

m_dot = [(launcher.st1.mp/launcher.st1.tb) (launcher.st2.mp/launcher.st2.tb) (launcher.st3.mp/launcher.st3.tb)];

tbb = [launcher.st1.tb launcher.st2.tb launcher.st3.tb];

L = [14.52 1.95 1.83];

Diam = 0.83;

thick = 0.001; %guess, must change. extract value from structural Masses

E = data.str.E; 
%%
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

LoD_time = Length./Diam;

Et_M = (E.*thick)./mass;

f1 = (9.87./(2.*pi)).*sqrt((pi.*Et_M)./(8.*(LoD_time.^3)));

Rho_eq = 5*(4*M0)/(pi*sum(L)*(((Diam+(2*thick))^2) - (Diam^2)));

how_many_modes = 2;

[xl,Xnx,~] = FFbeam(Diam,sum(L),thick,E,Rho_eq,how_many_modes);



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
