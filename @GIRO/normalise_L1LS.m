function [MZ_Adjustment, SamplesNormalised] = normalise_L1LS(OBJ_GIRO, lambdaL1Norm)

% This function compute a multiplicative MZ adjustment vector for each
% channel of LCMS sample using TIC.
%OBJ_GIRO.numChannels = 1;

BsplIntegMZ = [.0417 .4583 .4583 .0417];

% Initialise the normalisation images:
OBJ_GIRO.lambdaL1Norm = lambdaL1Norm;

MZ_Adjustment = ones(OBJ_GIRO.numSamples, OBJ_GIRO.sizeDyadicRT);

% Call the solver:
options.verbose = 0;

options.mode = 1;

options.order = -1;

options.progTol = .1;

options.maxIter = 20;

% Normalisation basis:            
minScale1D = 4;

maxScale1D = 6;

minScale2D  = 5; 

maxScale2D = 10;

% Cut-off of the 2D B-spline basis:
cutoffScale2D = 0;

% A vector specifying the channels used for normalisation:
indChannelNorm = 2;

numChannelsNorm = length(indChannelNorm);
    
MZWin_k = [1 370];

MZ_IND = [1 371];

% for j = 1 : numChannelsNorm
%     
%     MZWin_k = cell2mat(OBJ_GIRO.MZWin( indChannelNorm(j) ) );
%     
%     MZ_IND = [MZ_IND, MZ_IND(end)+MZWin_k(2)];
%     
% end

SamplesNormalised = zeros(OBJ_GIRO.sizeDyadicRT, MZ_IND(end) - MZ_IND(1), OBJ_GIRO.numSamples);

% Re-load in and concatonate the samples for normalisation:
offsetMZ = zeros(1,OBJ_GIRO.numSamples);
    
offsetRT = zeros(1,OBJ_GIRO.numSamples);
    
cd('~/Data/SWATH/hdf5/smo/MZ_-8/')
    
for j = 1 : OBJ_GIRO.numSamples
    
    datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(1)) '/%d/%d/%d/%d/cs'], -8, 5, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
      
    offset = double(h5readatt(cell2mat(OBJ_GIRO.fileName(j)),datasetname, 'Offset'));
        
    offsetMZ(j) = offset(1);
        
    offsetRT(j) = offset(2);
   
end

OBJ_GIRO.offsetRT = max(offsetRT);

offsetRT = OBJ_GIRO.offsetRT - offsetRT + 1;

OBJ_GIRO.offsetMZ = max(offsetMZ); 

offsetMZ = OBJ_GIRO.offsetMZ - offsetMZ + 1;

% Load in seamass coefs:
for j = 1 : OBJ_GIRO.numSamples
 
    MZ_counter = 1;

    for k = indChannelNorm
 
%    MZWin_k = cell2mat(OBJ_GIRO.MZWin( k ) );
    
    datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(k)) '/%d/%d/%d/%d/cs'], -8, 5, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
    
    SamplesNormalised(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2), MZ_IND(MZ_counter) : MZ_IND(MZ_counter+1)-1 , j) = double(h5read(cell2mat(OBJ_GIRO.fileName(j)),datasetname,[offsetMZ(j) offsetRT(j)],[MZWin_k(2)-MZWin_k(1) + 1, OBJ_GIRO.RTWinCTN_LevelN(2)-OBJ_GIRO.RTWinCTN_LevelN(1)+1])');
 
    end
        
%   MZ direction B-spline interpolation using a 4-tap FIR filter:
    SamplesNormalised(:,:,j) = filter(BsplIntegMZ, 1, SamplesNormalised(:,:,j), [], 2);
       
end
        
% Deform the samples:
indRT = repmat(1:OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.numSamples, 1); % Interpolating sites
    
DeformField = (OBJ_GIRO.BsplDictDeform * OBJ_GIRO.CP')';
    
indKnot = indRT - DeformField; % Knots
    
d = 3; % Order of the B-splines
     
% Initial deformation      
for j = 1 : OBJ_GIRO.numSamples

    [IntBspl, DIBspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d,0);

    SamplesNormalised(:,:,j) = IntBspl * SamplesNormalised(:,:,j);
    
end

N = 2;

numBspl = cell(numChannelsNorm , 2);

OBJ_GIRO.CN = cell(OBJ_GIRO.numSamples, numChannelsNorm);

OBJ_GIRO.CN(1,:) = {1};

% Vectorize the samples to normalise:

% Constructing the normalising basis:
[numBspl1D, numBspl2D, OBJ_GIRO.BsplDictNorm] = OBJ_GIRO.get_BsplDict_norm(minScale1D, maxScale1D, minScale2D, maxScale2D, MZ_IND(end)-MZ_IND(1), N);

numBspl = [{numBspl1D}, {numBspl2D}]; 
    
vecChannels = zeros(OBJ_GIRO.numSamples, size(SamplesNormalised,1) * size(SamplesNormalised,2));

for j = 1 : OBJ_GIRO.numSamples
     
    for k = 1 : numChannelsNorm
      
%        vecChannels(j , (OBJ_GIRO.sizeDyadicRT*(k-1) + OBJ_GIRO.sizeDyadicRT * sizeMZ_Channel * (i-2) + 1) : (OBJ_GIRO.sizeDyadicRT*k + OBJ_GIRO.sizeDyadicRT * sizeMZ_Channel * (i-2))) = sum(SamplesNormalised(:, (OBJ_GIRO.indChannel(i) + resMZ*(k-1)) : (OBJ_GIRO.indChannel(i) + resMZ*k - 1), j), 2);
        
        tmp = SamplesNormalised(:,:,j);

        tmp = tmp(:);
        
        vecChannels(j , :) = tmp';
        
    end
    
end
    
clear tmp


OBJ_GIRO.RefNorm = mean(log(Anscombe(vecChannels)), 1);

%    OBJ_GIRO.RefNorm = log(Anscombe(vecChannel_i(1,:)));
    
for j = 1 : OBJ_GIRO.numSamples

    OBJ_GIRO.TarNorm = vecChannels(j,:);
    
    CNVector = zeros(size(OBJ_GIRO.BsplDictNorm,2),1);
    
    OBJ_GIRO.indCN = ones(size(OBJ_GIRO.BsplDictNorm,2),1);
    
    CNVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_norm, CNVector, [lambdaL1Norm(1)*ones(sum(numBspl1D), 1); lambdaL1Norm(2)*ones(sum(numBspl2D), 1)], options);
    
    OBJ_GIRO.indCN = (CNVector ~= 0);
    
    CNVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_norm, CNVector, zeros(size(CNVector)), options);
    
    numBspl1D = cell2mat(numBspl(1));
    
    numBspl2D = cell2mat(numBspl(2)); 
        
        % Split the CNVector according to cutoffScale2D:
%        indDiff1D = 1 : sum(numBspl1D);
    indNorm2D = (sum(numBspl1D) + sum(numBspl2D(1:cutoffScale2D)) + 1) : (sum([numBspl1D; numBspl2D]));
    
%    OBJ_GIRO.CN(j, i) = {CNVector(indNorm2D)};
        
    OBJ_GIRO.CN(j) = {CNVector(indNorm2D)};
    
%        DiffExp = exp(OBJ_GIRO.BsplDictNorm(1: OBJ_GIRO.sizeDyadicRT, indDiff1D) * CNVector(indDiff1D));
    MZ_Adjustment(j,:) = exp(OBJ_GIRO.BsplDictNorm(1: OBJ_GIRO.sizeDyadicRT, indNorm2D) * CNVector(indNorm2D))'; 
        
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

MZ_Adjustment = MZ_Adjustment(:, OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2));

SamplesNormalised = SamplesNormalised(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2),:,:);

for j = 1 : OBJ_GIRO.numSamples
    
    SamplesNormalised(:,:,j) = SamplesNormalised(:,:,j) .* repmat(MZ_Adjustment(j,:)', 1, size(SamplesNormalised,2));
   
end


end

