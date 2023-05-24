function rx=ode_solved(RI,VI,MU,A,TimeCase,TOF,graph)

% ode_solved.m - Solution of the 2 body problem given initial conditions,
%                defined in a proper function (ode_2bproblem),using ode113
%
% PROTOTYPE: rx=ode_solved(RI,VI,MU,A,TimeCase,TOF)
% 
% INPUT:
%       RI[3]           Initial position vector
%       VI[3]           Initial velocity vector
%       MU[1]           Gravitational constant of the body
%       A[1]            Semi-major axis
%       TimeCase [1]    Logical variable to specify the period of integration
%                           0 : plot from initial position to the position reached
%                               after a given time equals to TOF
%                           1 : plot over the total period
%       TOF [1]         Time for a given transfer
%       graph [logical] Logical value for plotting if is it 'true' plot also the orbit 
%
% OUTPUT:
%       rx [...x3]      Matrix containing the three vectors that define the orbit in
%                       x,y,z coordinates
%
% CONTRIBUTORS:
% Schieppati Marco Simone
% 
% VERSIONS:
% 2019-11-18
%
 

T=2*pi*sqrt(A^3/MU);
if TimeCase==1
    tspan=0:70:T;
elseif TimeCase==0
    
    tspan=linspace(0,TOF,1000);
    
end
options=odeset( 'RelTol',1e-13, 'AbsTol', 1e-14);
y0x=[RI(1)
    RI(2)
    RI(3)
    VI(1)
    VI(2)
    VI(3)];
[~,rx]=ode113(@(t,y) ode_2bproblem(t,y,MU),tspan,y0x,options);


%% plotting 
if graph

plot3(rx(:,1),rx(:,2),rx(:,3),'LineWidth',1.5)

end