function [AMtime] = autoMatricesProtub(datcom,vars)

tic

n_hprot = length(vars.hprot);

%% datcom
for k = 1:2
    for n = 1:n_hprot
        
        datcom.xcg = vars.xcg(k);
        datcom.hprot = vars.hprot(n);
        createFor006(datcom);
        
        if k == 1
            % initialize joined empty matrix and state
            
            [CoeffsF, State] = datcomParser('full');
            State.hprot = vars.hprot;
            CoeffsE = struct();
            fn = fieldnames(CoeffsF);
            
            for f = 1:numel(fn)
                CoeffsE.(fn{f}) = zeros([size(CoeffsF.(fn{f})),n_hprot]);
            end
            
        else
            
            currentCoeffs = datcomParser();
            
            for f = 1:numel(fn)
                CoeffsE.(fn{f})(:,:,:,:,n) = currentCoeffs.(fn{f});
            end
            
            
        end
        
        
        
    end
end

%% Save joined empty .mat file
Coeffs = CoeffsE;
save('empty','State','Coeffs');

%%

delete('for003.dat', 'for004.dat', 'for005.dat', 'for006.dat', 'for009.dat',...
    'for010.dat', 'for011.dat', 'for012.dat')

AMtime = toc;