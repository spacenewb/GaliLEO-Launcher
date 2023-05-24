function createFor006(datcom)
%% recalling the variables
Mach = datcom.Mach;
Alpha = datcom.Alpha;
Beta = datcom.Beta;
Alt = datcom.Alt;

C1 = datcom.Chord1;
C2 = datcom.Chord2;
shape = datcom.shape;

xcg = datcom.xcg;
D = datcom.D;
r = D/2;
S = pi*D^2/4;
Lnose = datcom.Lnose;
OgType = datcom.OgType;
Lcenter = datcom.Lcenter;
Npanel = datcom.Npanel;
Phif = datcom.Phif;
Ler = datcom.Ler;
d = datcom.d;
zup_raw = datcom.zup_raw;
Lmaxu_raw = datcom.Lmaxu_raw;
H = datcom.Height;

Nm = length(Mach);
Na = length(Alpha);
Nb = length(Beta);
Nalt = length(Alt);
%% protuberance data
if isfield(datcom,'xprot') && datcom.hprot ~= 0 
    protub_flag = 1;
    % constants
    xprot = datcom.xprot;
    nloc = datcom.nloc;
    lprot = datcom.lprot;
    wprot = datcom.wprot;
    % variables
    hprot = datcom.hprot;
else
    protub_flag = 0;
end
%% datcom
Xle1 = Lcenter + Lnose - d - C1;
diffC = C1-C2;

switch shape
    case 'rect'
        Xle2 = Lcenter + Lnose - d - C2;
        
    case 'iso'
        Xle2 = Lcenter + Lnose - d - diffC/2;
        
    case 'parall'
        Xle2 = Lcenter + Lnose - d;
end

%%% Defining Fin Section
Zup = [zup_raw/C1, zup_raw/C2];
Lmaxu = [Lmaxu_raw/C1, Lmaxu_raw/C2];
Lflatu = [(C1 - 2*Lmaxu_raw)/C1, (C2 - 2*Lmaxu_raw)/C2];

%%% Creating for005.dat file from previous data
if ismac   % mac procedure
    fid  =  fopen(strcat(pwd, '/for005.dat'),'w+');
else
    fid  =  fopen(strcat(pwd, '\for005.dat'),'w+');
end

%%%%%%%%%%%% Flight Conditions
%%%% Beta
fprintf(fid, ' $FLTCON\r\n');
fprintf(fid, '  BETA = ');
fprintf(fid, '%.2f,\r\n', Beta(1));
%%%% Alt
fprintf(fid, '  ALT = ');
fprintf(fid, '%d', Nm);
fprintf(fid, '*');
fprintf(fid, '%d.,\r\n', Alt(1));
%%%% Nmach
fprintf(fid, '  NMACH = ');
fprintf(fid, '%d., \r\n', Nm);
%%%% Mach
fprintf(fid, '  MACH = ');
if Nm > 11
    for M = 1:11
        fprintf(fid, '%.2f,',Mach(M));
    end
    fprintf(fid,  ' \r\n MACH(12) = ');
    for M = 12:Nm
        fprintf(fid, '%.2f',Mach(M));
        fprintf(fid, ',');
    end
else
    for M = 1:Nm
        fprintf(fid, '%.2f,',Mach(M));
    end
end

fprintf(fid, '\r\n');
%%%% Nalpha
fprintf(fid, '  NALPHA = ');
fprintf(fid, '%d., \r\n', Na);
%%%% Alpha
fprintf(fid, '  ALPHA = ');
if Na > 9
    for a = 1:9
        fprintf(fid, '%.1f,', Alpha(a));
    end
    fprintf(fid,  ' \r\n  ALPHA(10) = ');
    for a = 10:Na
        fprintf(fid, '%.1f,', Alpha(a));
    end
else
    for a = 1:Na
        fprintf(fid, '%.1f,', Alpha(a));
    end
end

fprintf(fid, '$');

%%%%%%%%%%%% Reference Quantities
fprintf(fid, '\r\n $REFQ\r\n');
%%%% XCG
fprintf(fid, '  XCG = ');
fprintf(fid, '%.4f,\r\n', xcg);
%%%% SREF
fprintf(fid, '  SREF = ');
fprintf(fid, '%.5f,\r\n', S);
%%%% LREF
fprintf(fid, '  LREF = ');
fprintf(fid, '%.3f,\r\n', D);
%%%% LATREF
fprintf(fid, '  LATREF = ');
fprintf(fid, '%.3f,$', D);
%%%% TNOSE
%%%%%%%%%%%% Axisymmetric Body Geometry
fprintf(fid, '\r\n $AXIBOD\r\n');
%%%% TNOSE
fprintf(fid, '  TNOSE = %s, \r\n', OgType);
if strcmp(OgType,'POWER')
    fprintf(fid, '  POWER = %.3f, \r\n', datcom.NosePower);
end
%%%% LNOSE
fprintf(fid, '  LNOSE = ');
fprintf(fid, '%.3f, \r\n', Lnose);
%%%% DNOSE
fprintf(fid, '  DNOSE = ');
fprintf(fid, '%.3f, \r\n', D);
%%%% LCENTR
fprintf(fid, '  LCENTR = ');
fprintf(fid, '%.3f, \r\n', Lcenter);
%%%% DCENTR
fprintf(fid, '  DCENTR = ');
fprintf(fid, '%.3f, \r\n', D);
%%%% DEXIT
fprintf(fid, '  DEXIT = 0.,');
fprintf(fid, '\r\n');
%%%% BASE
fprintf(fid, '  BASE = .FALSE.,$');

% add protuberance if configuration is empty and if current
% %       protuberance area is different from zero
if protub_flag
    fprintf(fid, '\r\n $PROTUB \r\n');
    fprintf(fid, '  NPROT = 1., \r\n');
    fprintf(fid, '  PTYPE = BLOCK, \r\n');
    fprintf(fid, '  XPROT = ');
    fprintf(fid, '%.3f, \r\n', xprot);
    fprintf(fid, '  NLOC = ');
    fprintf(fid, '%.3f, \r\n', nloc);
    fprintf(fid, '  LPROT = ');
    fprintf(fid, '%.3f, \r\n', lprot);
    fprintf(fid, '  WPROT = ');
    fprintf(fid, '%.3f, \r\n', wprot);
    fprintf(fid, '  HPROT = ');
    fprintf(fid, '%.3f, \r\n', hprot);
    fprintf(fid, '  OPROT = 0.,$');
end

%%%%%%%%%%%% Finset
fprintf(fid, '\r\n $FINSET1 \r\n');
%%%% XLE
fprintf(fid, '  XLE = ');
fprintf(fid, '%.3f,', Xle1);
fprintf(fid, '%.3f, \r\n', Xle2);
%%%% NPANEL
fprintf(fid, '  NPANEL = ');
fprintf(fid, '%.1f,\r\n', Npanel);
%%%% PHIF
fprintf(fid, '  PHIF = ');
for P = 1:length(Phif)
    fprintf(fid, '%.1f', Phif(P));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
%%%% LER
fprintf(fid, '  LER = 2*');
fprintf(fid, '%.4f,\r\n', Ler);
%%%% SSPAN
fprintf(fid, '  SSPAN = ');
fprintf(fid, '%.3f,', r);
fprintf(fid, '%.3f,\r\n', r + H);
%%%% CHORD
fprintf(fid, '  CHORD = ');
fprintf(fid, '%.3f,', C1);
fprintf(fid, '%.3f,\r\n', C2);
%%%% SECTYP
fprintf(fid, '  SECTYP = HEX, \r\n');
%%%% ZUPPER
fprintf(fid, '  ZUPPER = ');
fprintf(fid, '%.4f,', Zup(1));
fprintf(fid, '%.4f,\r\n', Zup(2));
%%%% LMAXU
fprintf(fid, '  LMAXU = ');
fprintf(fid, '%.4f,', Lmaxu(1));
fprintf(fid, '%.4f,\r\n', Lmaxu(2));
%%%% LMAXU
fprintf(fid, '  LFLATU = ');
fprintf(fid, '%.4f,', Lflatu(1));
fprintf(fid, '%.4f,$ \r\n', Lflatu(2));

%%%%%%%%%%%% Options
fprintf(fid, 'DERIV RAD \r\n');
fprintf(fid, 'DIM M \r\n');
fprintf(fid, 'DAMP \r\n');
fprintf(fid, 'SAVE \r\n');
fprintf(fid, 'NEXT CASE \r\n');

%%%%%%%%%%%% Cases
for A = 1:Nalt
    for B = 1:Nb
        if not(A == 1 && B == 1)
            fprintf(fid,' $FLTCON \r\n');
            fprintf(fid,'  BETA = ');
            fprintf(fid, '%.1f, \r\n', Beta(B));
            fprintf(fid, '  ALT = ');
            fprintf(fid, '%d', Nm);
            fprintf(fid, '*');
            fprintf(fid, '%d.,$ \r\n', Alt(A));
            fprintf(fid, 'DERIV RAD\r\n');
            fprintf(fid, 'DIM M\r\n');
            fprintf(fid, 'DAMP\r\n');
            fprintf(fid, 'SAVE\r\n');
            fprintf(fid, 'NEXT CASE\r\n');
        end
    end
end
fclose(fid);

if ismac
    [~,~] = system('./datcom for005.dat' );
else
    [~,~] = system('datcom.exe for005.dat' );
end

value = 0;
while value == 0
    value = exist('for006.dat','file');
    pause(0.01);
end
