% function [RTAdjustment XCMS_Deformed] = BsplPyrRepDef_1D(obj)

% This function is a private member method for the class XCMS. It:
% 1. maintains the B-spline pyramid representations of XCMS data:
% 2. queries into the (1). Criterion (2). Deformation field 
%    (3). Normalising factors to iteratively deform the XCMS.
% 
%   [RTAdjustment XCMS_Deformed] = BsplPyrRep_1D(obj)
%
% the output: 
%   RTAdjustment is the retention time adjustment in the unit specified by obj.resoluteRT;
%   XCMS_Deformed is the deformed and normalised XCMS samples.

%% Initialising the iteration: 
%  Basic parameters:

% Going back to a very coarse level, so that there is no mis-alignment:
AppLowestRes = 30; % Approximate lowest resolution. In future this should be carefully considerred as a parameter
numResLevels = ceil(log2(obj.sizeRT / AppLowestRes));

sizeRT = zeros(1,numResLevels);

% The index and size of level j XCMS representation:
for i = numResLevels : (-1) : 1
    
Ind(i) = {[1 : 2^i : obj.sizeRT obj.sizeRT+1]};
    
sizeRT(i) = length(cell2mat(Ind(i)))-1;
    
end

% B-spline look-up-table:
CoarsestLengthSection = 2; % 8 pixel-wide B-splines
LengthBsplLUT = CoarsestLengthSection * 2^numResLevels;

u = (0 : LengthBsplLUT ) / LengthBsplLUT;

u = u(1:(end-1));

Bu3 = u.^3/6;
Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;
Bu0 = (1-u).^3/6;

LUT_Bu = [Bu3 Bu2 Bu1 Bu0];

% Gradient of B-splines:
dBu3 = u.^2/2;
dBu2 = -1.5*u.^2 + u + .5;
dBu1 = 1.5*u.^2 - 2*u;
dBu0 = -(1-u).^2/2;

LUT_dBu = [dBu3 dBu2 dBu1 dBu0];



% Interpolate the normalising field in the log domain to initialise the
% deformation in a finer scale:
% DiffScale = 3;
% 
% IniDeformXCMS = zeros(sizeRT(numResLevels-DiffScale), obj.sizeMZ, obj.numSamples);
% 
% ind = cell2mat(Ind(numResLevels-DiffScale));
% 
% for i = 1 : sizeRT(numResLevels-DiffScale)
%         
%     IniDeformXCMS(i,:,:) = sum(obj.XCMS_RAW(ind(i):ind(i+1)-1 , : , :), 1); % Initialize the images at level j
%     
% end
% 
% IterNFactor = zeros(size(IniDeformXCMS));
% 
% IterNFactor(1:2^DiffScale:end,:,:) = IniNFactor; 
% 
% IterNFactor(end,:,:) = IniNFactor(end,:,:);
% 
% for i = 1:obj.numSamples
%     
%     for j = 1:obj.sizeMZ
%     
%     IterNFactor(:,j,i) = interp1([1:2^DiffScale:sizeRT(numResLevels-DiffScale) sizeRT(numResLevels-DiffScale)], [IniNFactor(:,j,i); IniNFactor(end,j,i)], 1:sizeRT(numResLevels-DiffScale), 'linear');
%     
%     end
%     
% end
%     
% % Derive the multiplicative factor for the normalisation:
% IterNFactor = exp(IterNFactor);    

% Upsample and interpolate to make the NFactorIter applicable to a finer
% resolution:

% Normalisation: toward the geometirc mean of the samples: Initially
% applied to the images before registration to see if the criterion
% goes down.
   
% Here ends the Initial normalisation. 
% Iterative deformation / normalisation starts in the next section.

%% Iterative normalisation-deformation:

NFlag = 1;
    
RFlag = 1;
        
NRThresh = 1e+5;

CP = [];

gradCP = [];

L1NormCP = 0; % Initial L1 norm of control points

lambdaL1Deform = .1; 

sizeStep = .8;

for i = numResLevels : (-1) : 1
    
%% Normalisation:

% This part of algorithm normalise the Level i image: the normalisationo 
% field does not propagate, rather it is the L1-LS fit of the
% log-differences.
   
NormXCMS = zeros(sizeRT(i), obj.sizeMZ, obj.numSamples);
   
ind = cell2mat(Ind(i));

for j = 1 : sizeRT(i)
   
    NormXCMS(j,:,:) = sum(obj.XCMS_RAW(ind(j):ind(j+1)-1 , : , :), 1); % Initialize the images at level j
       
end
 
% Derive image normalisation: the difference of the log-image and the
% log-mean:

[ctnN gradCtn] = ctnVar2D(NormXCMS);

NormXCMS = NormXCMS .* NMASK; 

[ctnNNew gradCtn] = ctnVar2D(NormXCMS);

deltaCtn = ctnN - ctnNNew;
    
if deltaCtn > NRThresh
    
    NFlag = 1; % Flag of normalisation: to continue
 
else
    
    NFlag = 0;
    
end 



%% Registration:

indLevels = [0 cumsum(sizeRT(i:numResLevels)/(CoarsestLengthSection) + 3)];

MultiResBsplBasis = zeros(sizeRT(i), indLevels(end));


% Construct the Level i multi-resolution deformation field:
deltaCtn = 100000000; % Initialize the decrease of the criterion

ctnR = ctnN  + lambdaL1Deform*sum(abs(CP(:)));

dI_dg = zeros(size(gradCtn));

while RFlag 
    
    for j = 1 : obj.numSamples
        
        % Least square image representation:
    %     for i = 1 : obj.numSamples
%             
%         BsplInterpCoef(:,:,i) = BsplDict(1:sizeRT(j),1:sizeRT(j)) \ TMP_XCMS(:,:,i);    
%     
%         TMP_XCMS(:,:,i) = BsplDict(1:sizeRT(j),1:sizeRT(j)) * BsplInterpCoef(:,:,i);
%         
%     end
        
        % Query the image gradient for dI_dg:
        dI_dg(1:end-1,:,j) = diff(NormXCMS(:,:,j));
    
        dI_dg(end,:,j) = dI_dg(end-1,:,j);
        
        de_dg = sum(gradCtn(:,:,j) .* dI_dg(:,:,j) ./ (NormXCMS(:,:,j).^2), 2);
    
        % 
    
    % Chain rule to combine the three: integrate across the deformation 
    % field for dg_dC:
        gradCP(j,:) = de_dg' * MultiResBsplBasis;
    
    end
    
    % Shrinkage: adding L1 regularisation:
    gradCP = gradCP + lambdaL1*abs(CP);
        
    % Update the control points:
    CP = CP - sizeStep * gradCP;
    
    % Generating the deformation field:
    DeformField = (MultiResBsplBasis * CP')';
    
    % Normalising the deformation field: set the maximum displacement:
    DeformField = DeformField * (2 / max(DeformField(:)));
     
    % Volume-preserving interpolation:
    DeformXCMS = Interp1D_MassDist(NormXCMS, DeformField);    
        
    % Criterion re-evaluation:
    [ctnRNew gradCtnR] = ctnVar2D(DeformXCMS);
    
    ctnRNew = ctnRNew + lambdaL1Deform*sum(abs(CP(:)));
    
    deltaCtn = ctnR - ctnRNew; 
    
    if deltaCtn > NRThresh
        
        % Continue deforming
        RFlag = 1; 
        
    else
        
        %Stop deforming: 
        RFlag = 0;
        
        % reverse the previous control point update
        CP = CP + sizeStep * gradCP;
    
        % Generating the deformation field:
        DeformField = (MultiResBsplBasis * CP')'; 

    
    end

    
end
    
end

%     
% tmp = zeros(size(Iter_XCMS));
% 
% sizeRT_NFactorCum = size(NFactorCum);
%         
% tmp(1:2:2*sizeRT_NFactorCum,:,:) = NFactorCum;
%         
% tmp(2:2:2*sizeRT_NFactorCum,:,:) = NFactorCum;
%         
% tmp(end,:,:) = NFactorCum(end,:,:);
%         
% NFactorCum = tmp;


   
%%
% Image representation organised in data array
%
% for j = 1 : numResLevels
%     
%     sizeRT(j) = size(TMP_XCMS,1);
%     
%     for i = 1 : obj.numSamples
%             
%         BsplInterpCoef(:,:,i) = BsplDict(1:sizeRT(j),1:sizeRT(j)) \ TMP_XCMS(:,:,i);    
%     
%         TMP_XCMS(:,:,i) = BsplDict(1:sizeRT(j),1:sizeRT(j)) * BsplInterpCoef(:,:,i);
%         
%     end
% 
%     % Image:
%     XCMS_BsplPyrRep(j) = {TMP_XCMS};
%         
%     % Coefficients
%     BsplPyrRepCoef(j) = {BsplInterpCoef};
% 
%     % Downsampling
%     TMP_XCMS = TMP_XCMS(1:2:end,:,:);
%     
%     BsplInterpCoef = BsplInterpCoef(1:2:end,:,:);
%     
% end
% 
% clear TMP_XCMS BsplInterCoef


%% Deforming / Normalising the samples iteratively:
% 
% for j = numResLevels : (-1) : 1
%
%     tmp_xcms = Iter_XCMS;
%         
%     deltaCtn = 100000000; % Initialize the decrease of the criterion
%     
%     [ctnN gradCtn] = ctnVar2D(Iter_XCMS);
%     
%     ctnR = ctnN  + lambdaL1*sum(abs(CP(:)));
%     
%     % Applying normalising factor after the second iteration
%     if j == numResLevels 
%       
%     mu = .25;
%     
%     dampMu = 1;%.85;
%     
% while (NFlag && RFlag) % Stopping criterion
%     
%   % %    NFactorTmp = NFactorCum .* NFactorIter;
%      
%     NormalisedXCMS = Iter_XCMS .* NFactorIter; 
%     
%     [ctnNNew gradCtn] = ctnVar2D(NormalisedXCMS);
%     
%     deltaCtn = ctnN - ctnNNew;
%     
%     if deltaCtn > NRThresh
%         
%         Iter_XCMS = NormalisedXCMS;
%         
%         NFactorCum = NFactorCum .* NFactorIter;
%                 
%         tmp = tmp_xcms .* NFactorCum;
%     
%         sum(sum(sum(NormalisedXCMS - tmp)))
%         
%         ctnN = ctnNNew;
%         
%         NFlag = 1; % Flag of normalisation: to continue
%     
%         mu = mu * dampMu;
%         
%     else
%         
%         NFlag = 0;
%         
%     end 
%     
%     % Multiresolution deformation field associated with level j:     
%     MultiResBsplDeformBasis = zeros(sizeRT(j), sum(sizeRT(j:numResLevels)));
%     indLevels = [0 cumsum(sizeRT(j:numResLevels))];
%     
%     for i = j : numResLevels
%     
%     % Downsample the look-up-table:
%     downFactor = 2^(i-j);
%     Bu = LUT_Bu(1:(LengthB/downFactor):end); % Length of Bu is a multiple of 4
%     
%     Bu = Bu / max(Bu); % L_inf normalise
%     
%     LevelJDeformField = convmtx(Bu, sizeRT(j));
% 
%     LevelJDeformField = LevelJDeformField(1:downFactor:end, (length(Bu)/2+1) : (length(Bu)/2 + sizeRT(j))); 
%     
%     MultiResBsplDeformBasis(:,indLevels(i-j+1)+1:indLevels(i-j+2)) = LevelJDeformField';
%         
%     end
% 
%     MultiResBsplDeformBasis = sparse(MultiResBsplDeformBasis); % Deformation field dictionary
%      
%     CP = [zeros(obj.numSamples, sizeRT(j)), CP]; % Control Point
%                     
%     gradCP = [zeros(obj.numSamples, sizeRT(j)), gradCP]; % Gradient for control point update
% 
%      
%     BsplPyrRepCoef = zeros(sizeRT(j), obj.sizeMZ, obj.numSamples);
%     
%     for i = 1 : obj.numSamples
%         
%         tic;
%             
%         for m = 1 : obj.sizeMZ
%                     
%         BsplPyrRepCoef(:,m,i) = lsqnonneg(BsplDict(1:sizeRT(j),1:sizeRT(j)) , Iter_XCMS(:,m,i));    
%     
%         end
%          
%         toc
% 
%         sum(sum(abs(Iter_XCMS(:,:,i)-BsplDict(1:sizeRT(j),1:sizeRT(j))*BsplPyrRepCoef(:,:,i))))
%         
%     end
%     
%     while RFlag
%     
%         TIC = zeros(obj.numSamples,sizeRT(j));
%         
%         % Chain rule:
%         for i = 1 : obj.numSamples 
%         
%         TIC(i , :) = sum(Iter_XCMS(:,:,i),2)';    
%             
%         gradImg = (gradBsplDict(1:sizeRT(j), 1:sizeRT(j)) * BsplPyrRepCoef(:,:,i));  
%                     
%         NormI = (Iter_XCMS(:,:,i).^2);
%         
%         NormI(NormI < 1) = 1;
%         
%         de_dI = sum( (gradCtn(:,:,i) .* gradImg) ./ NormI , 2)';
%                 
% %        de_dI = de_dI ./ max(de_dI);
%          
%         gradCP(i,:) = de_dI * MultiResBsplDeformBasis;
%         
%         gradCP(i,:) = gradCP(i,:) ./ max(abs(gradCP(i,:)));
%         
%         end 
%         
%         % Shrinkage:
%         gradCP = gradCP + lambdaL1*abs(CP);
%         
%         % Update the control points:
%         CP = sizeStep * gradCP;
%         
%         DeformField = (MultiResBsplDeformBasis * CP')';
%         
%         % Interpolating the XCMS to get the deformed XCMS:
%         XCMS_Deformed = zeros(size(Iter_XCMS));
%         
%         for i = 1 : obj.numSamples
%          
%             BsplDictDeformed = zeros(sizeRT(j),sizeRT(j));
%             
%             BsplDictDeformed([1:4 sizeRT(j)-4:sizeRT(j)],:) = BsplDict([1:4 sizeRT(j)-4:sizeRT(j)],1:sizeRT(j)); 
%             
%             for m = 5 : sizeRT(j)-5
%                 
%                 % New index for the knots:
%                 indNew = ((m-2) : (m+2)) + DeformField(i,((m-2) : (m+2)));
%             
%                 % Interpolate the integers in between 1st segment:
%                 if (floor(indNew(2)) - floor(indNew(1))) >= 1
%                 % Look-up:
%                 indInterp1 = ceil(indNew(1)) : floor(indNew(2));
%                 BsplDictDeformed(indInterp1,m) = LUT_Bu(1 + round((indInterp1 - indNew(1)) * LengthB / (indNew(2)-indNew(1))));
%                 end
%                 
%                 if (floor(indNew(3)) - floor(indNew(2))) >= 1
%                 % Interpolate the integers in between 2nd segment:
%                 indInterp2 = ceil(indNew(2)) : floor(indNew(3));
%                 % Look-up:
%                 BsplDictDeformed(indInterp2,m) = LUT_Bu(LengthB + round((indInterp2 - indNew(2)) * LengthB / (indNew(3)-indNew(2))));
%                 end
%                 
%                 if (floor(indNew(4)) - floor(indNew(3))) >= 1
%                 % Interpolate the integers in between 3rd segment:
%                 indInterp3 = ceil(indNew(3)) : floor(indNew(4));
%                 % Look-up:
%                 BsplDictDeformed(indInterp3,m) = LUT_Bu( 2*LengthB + round((indInterp3 - indNew(3)) * LengthB / (indNew(4)-indNew(3))));
%                 end
%                 
%                 if (floor(indNew(5)) - floor(indNew(4))) >= 1
%                 % Interpolate the integers in between 4th segment:
%                 indInterp4 = ceil(indNew(4)) : floor(indNew(5));
%                 % Look-up:
%                 BsplDictDeformed(indInterp4,m) = LUT_Bu( 3*LengthB + round((indInterp4 - indNew(4)) * LengthB / (indNew(5)-indNew(4))));
%                 end
%                 
%                 % Normalise:
%                 BsplDictDeformed(:,m) = BsplDictDeformed(:,m) * normBu / sum(BsplDictDeformed(:,m));
%                 
%             end
%             
%             % Derive the Deformed image:
%             XCMS_Deformed(:,:,i) = BsplDictDeformed * BsplPyrRepCoef(:,:,i);
%             
%         end
%              
%         [ctnNew gradCtn] = ctnVar2D(XCMS_Deformed);
%         
%         deltaCtn = ctn - ctnNew
%    
%         if deltaCtn > 0
%             
%             IterXCMS = XCMS_Deformed;
%             
%             RFlag = 1;
%             
%         else
%             
%             RFlag = 0;
%             
%         end
%         
%     end
%         
% end
%     
% 
%     
% end
% 
% 
%  
% RTAdjustment = zeros(obj.numSamples, obj.sizeRT);
% 
% XCMS_Deformed = zeros(size(XCMS));

%% 















