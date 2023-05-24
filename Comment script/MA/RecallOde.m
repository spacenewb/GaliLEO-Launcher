function [all_steps] = RecallOde(fun, T, Y, varargin)
%{

RECALLODEFCN - This function allows to compute some parameters used
inside the ODE integrations

INPUTS:
            - fun, ode function used;
            - T, integration time vector;
            - Y, state matrix.

OUTPUTS:
            - all_steps, structure which contains all the parameters needed.

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept
email: adriano.filippo.inno@skywarder.eu
Release date: 16/11/2018

%}

NT = length(T);
fun_info = functions(fun);

for i = 1:NT
    
    [~,single_step] = fun(T(i),Y(i,:),varargin{:});
    
    all_steps.air.rho(i) = single_step.air.rho;
    all_steps.air.P(i) = single_step.air.P;
    
    all_steps.accelerations.body_acc(i) = single_step.accelerations.body_acc;
    
    all_steps.interp.M(i) = single_step.interp.M;
    all_steps.interp.alpha(i) = single_step.interp.alpha;
%    all_steps.interp.vel(i) = single_step.interp.vel;
    
    all_steps.forces.AeroDyn_Forces(i, 1:3) = single_step.forces.AeroDyn_Forces;
%    all_steps.forces.T(i) = single_step.forces.T;
    
    
     all_steps.coeff.CD(i) = single_step.coeff.CD;
    % all_steps.coeff.XCP(i) = single_step.coeff.XCP;
    %
    %
end