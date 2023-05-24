function [Coeffs, State] = datcomParser(varargin)
%{ 

DATCOMPARSER - function to parse the 'for006.dat' file and to generate the
               matrices of the aerodynamic coefficients.

INPUTS:
            - varargin: it can be either a string that contains the name of
                        .mat file that is saved in the Current Folder or it
                        can be an empty field and the .mat file is not
                        saved. (Read NOTES for more informations.)

OUTPUTS: 
            - Coeffs:   struct in which the fields are the matrices 
                        of the aerodynamic coefficients. (e.g. Coeffs.CA)
            - State:    struct in which the fields are the arrays of the
                        aerodynamic states (Alpha, Mach, Beta, Altitude)
                        used to compute the aerodynamic coefficients.

NOTES: If the function is used with the syntax 
'[Coeffs, State] = datcomParser()' the file .mat is not saved in the 
Current Folder. If the function is used with the syntax 
'[Coeffs, State] = datcomParser('name')' a file called 'name.mat' is saved
in the Current Folder. This file contains the two structs Coeffs and
State that can be loaded into the workspace in any future moment using the
command MATLAB function 'load'. Note that the file 'for006.dat' that has to
be parsed must be in the Current Folder.

%}

if not(isempty(varargin))
    mat_name = varargin{1};
    savemat = true;
else
    savemat = false;
end

linestring = fileread('for006.dat');

%% blocksplit
pattern = '\*+ \<FLIGHT CONDITIONS AND REFERENCE QUANTITIES \*+';
blocks = regexp(linestring,pattern,'split');
block1 = blocks(1);
blocks = blocks(2:end);

%% error check
error_check = regexp(block1,'(** ERROR **)','tokens');
if not(isempty(error_check{1}))
    error('Attention: Error in the ''for006.dat'' file.');
end

%% get_coeffs_name
pattern =  ' *\<ALPHA\> *([\w -/]*)';
names = cell(26,1);
index = 1;

for i = 1:4
    block = blocks{i};
    token = regexp(block,pattern,'tokens');
    
    % convert cells inside token into strings
    for k = 1:length(token)
        token{k} = char(token{k});
    end
    
    dummy = strjoin(token);
    dummy = split(dummy);
    
    % replaceBadChars
    correct = dummy;
    pattern1 = '[\./-]';
    for j = 1:length(dummy)
        name = dummy{j};
        a = regexprep(name(1:end-1),pattern1,'_');
        b = regexprep(name(end),pattern1,'');
        correct{j} = [a,b];
    end
    
    names(index:index+length(correct)-1) = correct;
    index = index+length(correct);
end

%% get_data
pattern1 = ' *\<MACH NO\> *= * (-*[\d]+.[\d]*)';
pattern2 = ' *\<ALTITUDE\> *= * (-*[\d]+.[\d]*)';
pattern3 = ' *\<SIDESLIP\> *= * (-*[\d]+.[\d]*)';

M=cell(1,length(blocks)/4);
A=cell(1,length(blocks)/4);
B=cell(1,length(blocks)/4);

for i = 1:length(blocks)/4
    block = blocks{(i-1)*4+1};
    mach = regexp(block,pattern1,'tokens');
    M{i} = char(mach{1});
    sslip = regexp(block,pattern3,'tokens');
    B{i} = char(sslip{1});
    alt = regexp(block,pattern2,'tokens');
    A{i} = char(alt{1});
end

M = str2num(strjoin(M));
A = str2num(strjoin(A));
B = str2num(strjoin(B));
%% get_alpha
pattern = '^[-\d](?=[\d\.])';
pattern2 = '\n\t?';

block = blocks{2};
lines = regexp(block,pattern2,'split');
index = 0;
new_1 = cell(200,1);

for j = 1:length(lines)
    line = lines{j};
    line = strip(line);
    
    if regexp(line,pattern,'once')
        index = index + 1;
        new_1{index} = str2double(split(line));
    end
    
end

alpha = zeros(1,index);

for j = 1:index
    row = new_1{j};
    alpha(j) = row(1);
end
%% get_rawdata
raw_data = cell(1,length(blocks));

for i = 1:length(blocks)
    block = blocks{i};
    lines = regexp(block,pattern2,'split');
    index = 0;
    new_1 = cell(200,1);
    
    for j = 1:length(lines)
        line = lines{j};
        line = strip(line);
        
        if regexp(line,pattern,'once')
            index = index + 1;
            new_1{index} = transpose(str2num(line));
            
        end
        
    end
    
    new_1_1 = cell(1,index);
    
    for j = 1:index
        row = new_1{j};
        new_1_1{j} = row(2:end);
    end
    
    new_1 = new_1_1;
    l = length(new_1{1});
    
    for j = 1:length(new_1_1)
        row = new_1_1{j};
        if length(row)~=l
            new_1{j-length(new_1_1)/2} = [ new_1{j-length(new_1_1)/2}; new_1_1{j}];
            new_1{j} = {};
        else
            index = j;
        end
    end
    
    new_1 = transpose(cell2mat(new_1(1:index)));
    
    raw_data{i} = new_1;
    
end


raw_data = cell2mat(raw_data);

%% savemat
realM = [M(1), NaN(1,200)];
realA = [A(1), NaN(1,200)];
realB = [B(1), NaN(1,200)];

iM = 1;
iA = 1;
iB = 1;

for i = 2:length(M)
    if not(any(realM == M(i)))
        iM = iM + 1;
        realM(iM) = M(i);
    end
    if not(any(realA == A(i)))
        iA = iA + 1;
        realA(iA) = A(i);
    end
    if not(any(realB == B(i)))
        iB = iB + 1;
        realB(iB) = B(i);
    end
end

realM = realM(1:iM);
realA = realA(1:iA);
realB = realB(1:iB);

for j = 1:length(names)
    Coeffs.(names{j}) = zeros(length(alpha),iM,iB,iA);
end

for i = 1:length(blocks)/4
    index = i;
    iA = realA==A(index);
    iB = realB==B(index);
    iM = realM==M(index);
    
    for j = 1:length(names)
        Coeffs.(names{j})(:,iM,iB,iA) = raw_data(:,length(names)*(i-1)+j);
    end
    
    
end

State.Machs = realM;
State.Alphas = alpha;
State.Betas = realB;
State.Altitudes = realA;

if savemat
    save(mat_name,'State','Coeffs');
end

end

