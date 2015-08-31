function [MZ_Adjustment, SamplesNormalised] = normalise_L1LS(OBJ_GIRO, lambdaL1Norm)

% This function compute a multiplicative MZ adjustment vector for each
% channel of LCMS sample.
%OBJ_GIRO.numChannels = 1;

% Initialise the normalisation images:
OBJ_GIRO.lambdaL1Norm = lambdaL1Norm;

SamplesNormalised = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples); 

SamplesNormalised(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2), :, :) = OBJ_GIRO.SamplesDeformed;

MZ_Adjustment = ones(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.numChannels, OBJ_GIRO.numSamples);

% Call the solver:
options.verbose = 0;

options.mode = 1;

options.order = -1;

options.progTol = .1;

options.maxIter = 80;

% Normalisation basis:            
minScale1D = 3;

maxScale_1D = 3;

minScale2D  = 4;

maxScale2D = 10;

% Cut-off of the 2D B-spline basis:
cutoffScale2D = 0;

N = 2;

numBspl = cell(OBJ_GIRO.numChannels , 2);

OBJ_GIRO.CN = cell(OBJ_GIRO.numSamples, OBJ_GIRO.numChannels);

OBJ_GIRO.CN(1,:) = {1};

for i = 1 : OBJ_GIRO.numChannels
    
    sizeMZ_Channel = OBJ_GIRO.indChannel(i+1) - OBJ_GIRO.indChannel(i);
    
    % Constructing the normalising basis:
    [numBspl1D, numBspl2D, OBJ_GIRO.BsplDictNorm] = OBJ_GIRO.get_BsplDict_norm(minScale1D, maxScale_1D, minScale2D, maxScale2D, sizeMZ_Channel, N);

    numBspl(i, 1) = {numBspl1D};
    
    numBspl(i, 2) = {numBspl2D}; 
    
    % Load in a channel of data and vectorize it:
    vecChannel_i = zeros(OBJ_GIRO.numSamples, OBJ_GIRO.sizeDyadicRT * sizeMZ_Channel);
    
    for j = 1 : OBJ_GIRO.numSamples
        
        tmp = SamplesNormalised(:, OBJ_GIRO.indChannel(i) : (OBJ_GIRO.indChannel(i+1) - 1), j);
        
        vecChannel_i(j,:) = tmp(:)';
        
    end
    
    clear tmp

    OBJ_GIRO.RefNorm = mean(log(Anscombe(vecChannel_i)), 1);

%    OBJ_GIRO.RefNorm = log(Anscombe(vecChannel_i(1,:)));
    
    for j = 1 : OBJ_GIRO.numSamples
    
        OBJ_GIRO.TarNorm = vecChannel_i(j,:);
    
        CNVector = zeros(size(OBJ_GIRO.BsplDictNorm,2),1);
        
        OBJ_GIRO.indCN = ones(size(OBJ_GIRO.BsplDictNorm,2),1);
        
        CNVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_norm, CNVector, lambdaL1Norm*ones(size(CNVector)), options);
     
        OBJ_GIRO.indCN = (CNVector ~= 0);
        
        CNVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_norm, CNVector, zeros(size(CNVector)), options);
        
        numBspl1D = cell2mat(numBspl(i, 1));
    
        numBspl2D = cell2mat(numBspl(i, 2)); 
        
        % Split the CNVector according to cutoffScale2D:
%        indDiff1D = 1 : sum(numBspl1D);
        
        indNorm2D = (sum(numBspl1D) + sum(numBspl2D(1:cutoffScale2D)) + 1) : (sum([numBspl1D; numBspl2D]));
        
        OBJ_GIRO.CN(j, i) = {CNVector(indNorm2D)};
        
%        DiffExp = exp(OBJ_GIRO.BsplDictNorm(1: OBJ_GIRO.sizeDyadicRT, indDiff1D) * CNVector(indDiff1D));
        
        MZ_Adjustment(:, i, j) = exp(OBJ_GIRO.BsplDictNorm(1: OBJ_GIRO.sizeDyadicRT, indNorm2D) * CNVector(indNorm2D)); 
        
%         figure; plot(OBJ_GIRO.RefNorm,'b');
%         
%         hold on;
%         
%         plot(log(Anscombe(OBJ_GIRO.TarNorm)),'r')
%         
%         plot(log(Anscombe(OBJ_GIRO.TarNorm .* repmat(MZ_Adjustment(:, i, j), length(OBJ_GIRO.TarNorm) / length(MZ_Adjustment(:,i,j)), 1)')),'g');
%         
%         tmp = exp(OBJ_GIRO.BsplDictNorm * CNVector); 
%         
%         plot(log(Anscombe(OBJ_GIRO.TarNorm .* tmp')), 'y');
        
    end
        
end

MZ_Adjustment = MZ_Adjustment(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2),:,:);

SamplesNormalised = SamplesNormalised(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2),:,:);

end

