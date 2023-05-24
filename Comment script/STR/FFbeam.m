function [xl,Xnx,fn] =  FFbeam(Diam,Length,thickness,E,Rho,Mn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Modal frequencies and the Modal shapes for a Free-Free Beam model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants Declaration

R=Diam/2;
L=Length;
Ix=pi*thickness*R^3;  % cylinder with ring cross section
A=pi*(((R+thickness)^2) - (R^2));
 
HMMS=Mn; % The number of modes and mode shapes to be computed


%% Calculate the Beta (characteristic shape) values for each mode

Nm=3*HMMS;
    jj=1;
    
    while jj<=Nm
        betaNL(jj)=fzero(@(betaNL)cos(betaNL)*cosh(betaNL)-1,jj+3);
        jj=jj+3;
    end
    
    index=(betaNL~=0);
    betaNLall=(betaNL(index))'; % Beta values
    betaN=(betaNLall/L)';   % Normalised Beta values
    
%% Compute the natural frequencies (fn) of the launcher

k=1;
wn=zeros(1,length(betaN));  % Angular natural frequency
fn=ones(1,length(wn));  % Natural frequency in Hz

while k<=length(betaN)
    wn(k)=betaN(k)^2*sqrt((E*Ix)/(Rho*A)); 
    fn(k)=wn(k)/(2*pi); 
    %fprintf('Mode shape # %2f corresponds to nat. freq (fn): %3.3f\n', k, fn(k) )
    k=k+1;
end


%% Declare and compute the Sigma Function for each Beta

x=linspace(0, L, 720); % create 720 elements in the entire beam (arbitrary value 720)
sigmaN=zeros(1,HMMS);

for ii=1:HMMS
    sigmaN(ii)=(cosh(betaNLall(ii))-cos(betaNLall(ii)))/(sinh(betaNLall(ii))-sin(betaNLall(ii)));
end

xl=x./L; % Normalised location along the beam for an element
Tc='(cosh(betaN(ii).*x(jj))+cos(betaN(ii).*x(jj)))-sigmaN(ii).*(sinh(betaN(ii).*x(jj))+sin(betaN(ii)*x(jj)))';  % Complete shape function

Xnx=zeros(length(betaN),length(x));

%% Solving the shape functions for each Beta
for ii=1:length(betaN)
    for jj=1:length(x)
        Xnx(ii,jj)=eval(Tc); % Displacement of each element
    end
end

XnxMAX=max(abs(Xnx(1,1:end))); % Maximum displacement
Xnx=Xnx./XnxMAX; % Normalisation of the displaceents
Xnx=(Xnx./L).*((xl.*1.5)+1); % Mapping of the displaceents to the beam model

end
