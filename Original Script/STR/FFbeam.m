function [xl,Xnx,fn] =  FFbeam(Diam,Length,thickness,E,Rho,Mn)

R=Diam/2;
L=Length;
Ix=pi*thickness*R^3; %% cylinder with ring cross section
A=pi*(((R+thickness)^2) - (R^2));
 
HMMS=Mn; % The number of modes and mode shapes to be computed

if HMMS>=7
    disp('   ')
    warning('NOTE: Up to 6 mode shapes (plots) are displayed via the script. Yet, using evaluated data (Xnx) of the script, more mode shapes can be plotted');
    disp('   ')
end

    Nm=3*HMMS;
    jj=1;
    
    while jj<=Nm
        betaNL(jj)=fzero(@(betaNL)cos(betaNL)*cosh(betaNL)-1,jj+3);
        jj=jj+3;
    end
    
    index=(betaNL~=0);
    betaNLall=(betaNL(index))';
    betaN=(betaNLall/L)';
    
k=1;
wn=zeros(1,length(betaN));
fn=ones(1,length(wn));

while k<=length(betaN)
    wn(k)=betaN(k)^2*sqrt((E*Ix)/(Rho*A));
    fn(k)=wn(k)/(2*pi);
    %fprintf('Mode shape # %2f corresponds to nat. freq (fn): %3.3f\n', k, fn(k) )
    k=k+1;
end

x=linspace(0, L, 720);
sigmaN=zeros(1,HMMS);

for ii=1:HMMS
    sigmaN(ii)=(cosh(betaNLall(ii))-cos(betaNLall(ii)))/(sinh(betaNLall(ii))-sin(betaNLall(ii)));
end

xl=x./L;
Tc='(cosh(betaN(ii).*x(jj))+cos(betaN(ii).*x(jj)))-sigmaN(ii).*(sinh(betaN(ii).*x(jj))+sin(betaN(ii)*x(jj)))';
Xnx=zeros(length(betaN),length(x));

for ii=1:length(betaN)
    for jj=1:length(x)
        Xnx(ii,jj)=eval(Tc);
    end
end

XnxMAX=max(abs(Xnx(1,1:end)));
Xnx=Xnx./XnxMAX;
Xnx=(Xnx./L).*((xl.*1.5)+1);

end
