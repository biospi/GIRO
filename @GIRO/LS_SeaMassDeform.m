function RT_Adjustment = LS_SeaMassDeform(obj)

% This is an experimental version of the deformation algorithm, in the
% sense:
% 1. Data should be read and passed to the algorithm through another
% function, but in this experiment the data read algorithm is within the
% deformation algorithm. 
mz_resolution = -6;

shrinkage = -3;

tolerance = -6;
              
sizeNorm = 4;%[4 4 8 8 16 16 32 32]; % Size between two knots of B-spline basis in normalisation

BsplIntegMZ = [.0417 .4583 .4583 .0417];
        
% Iterative normalisation-deformation: 
% Normalising flag: 1 if a new normalisation reduces the variation more 
% than a threshold NThresh, which lead to new deformation; 0 otherwise and
% stop the normalisation, going to the next level:
% NFlag = 1; 
% NThresh = 400;

lambdaL1Deform = 0; 

ResBegin = 3;

ResEnd = 6;

%numResLevels = ResEnd - ResBegin + 1;

% Retention window used for computing the criterion, meaningful retention
% window at each level is derived from this interval:

for i = ResBegin : ResEnd
  
    offsetMZ = zeros(1,obj.numSamples);
    
    offsetRT = zeros(1,obj.numSamples);
    
    % Load in the SeaMass coefs:    
    for j = 1 : obj.numSamples
            
        datasetname = sprintf('/%d/%d/%d/%d/0/Coefs', mz_resolution, i, shrinkage, tolerance);
    
        % Offset in RT dimension is ignored because it will be corrected by deformation: 
        offset = double(h5readatt(obj.fileName(j,:),datasetname, 'Offset'));
        
        offsetMZ(j) = offset(1);
        
        offsetRT(j) = offset(2);
        
    end 
        
    RTWin_i = round(obj.RTWin / 2^(ResEnd - i)) + 1;
    
    sizeRT_i = RTWin_i(2) - RTWin_i(1) + 1;
    
    sizeDyadicRT = pow2(ceil(log2(sizeRT_i)));    
        
    obj.RTWinCTN_LevelN(1) = round((sizeDyadicRT - sizeRT_i) / 2);
    
    obj.RTWinCTN_LevelN(2) = obj.RTWinCTN_LevelN(1) + sizeRT_i - 1;
    
    RelativeOffsetMZ = offsetMZ - min(offsetMZ);
    
    % Allocate a dyadically sized RT for the samples:
    SeaMassCoef = zeros(sizeDyadicRT, obj.sizeMZ, obj.numSamples);
    
    for j = 1 : obj.numSamples
         
        datasetname = sprintf('/%d/%d/%d/%d/0/Coefs', mz_resolution, i, shrinkage, tolerance);
        
        SeaMassCoef(obj.RTWinCTN_LevelN(1):obj.RTWinCTN_LevelN(2),:,j) = double(h5read(obj.fileName(j,:),datasetname,[1+RelativeOffsetMZ(j) offsetRT(j)],[obj.sizeMZ sizeRT_i])');
                                       
        % MZ direction B-spline interpolation using a 4-tap FIR filter:
        SeaMassCoef(:,:,j) = filter(BsplIntegMZ, 1, SeaMassCoef(:,:,j), [], 2);
            
    end
        
    obj.SeaMassCoef = SeaMassCoef;
    
    %% Establish the level i multiresolution deformation dictionary:
    [obj.BsplDictDeform normBuL2 numCoefs] = obj.BsplDeformDict(sizeDyadicRT, i - ResBegin + 1);
              
    % and the corresponding control points:
    % Attach the new control points to the control points from 
    % previous levels: 
    if i - ResBegin >= 1      
    
     indLevels = [0 cumsum(numCoefs(2:end))];   
                 
    for j = 1 : (i-ResBegin) 
        
        obj.CP(:, indLevels(j) + 1 : indLevels(j+1)) = obj.CP(:, indLevels(j) + 1 : indLevels(j+1)) * 2 * normBuL2(j+1) / normBuL2(j);
        
    end
    
    end 
     
    % Taken as initial of the next level:
    obj.CP = [zeros(obj.numSamples, size(obj.BsplDictDeform,2) - size(obj.CP,2)), obj.CP]; 

    %% Construct the integrated 3rd order B-spline basis matrix and its 
    %  derivative matrix: the 3rd order B-spline basis matrix
    indRT = repmat(1:sizeDyadicRT, obj.numSamples, 1); % Interpolating sites
    
    DeformField = (obj.BsplDictDeform * obj.CP')';
    
    indKnot = indRT - DeformField; % Knots
    
    d = 3; % Order of the B-splines
     
    XCMS_Deform = zeros(sizeDyadicRT, obj.sizeMZ, obj.numSamples);
    
    % Initial deformation      
    for j = 1 : obj.numSamples
    
        [IntBspl Bspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d,0);

        XCMS_Deform(:,:,j) = IntBspl * SeaMassCoef(:,:,j);
         
    end
    
     
%    while NFlag ~= 0    
    %% Normalising:
    obj.NMASK = obj.NormaliseLS(XCMS_Deform, sizeNorm*2^(i-ResBegin+1)); %sizeNorm*2^(i+1)
 
 
        
    %% Deforming:
        
%   % Call the solver:
    options.verbose = 0;
   
    options.mode = 1; 
  
    options.order = 1;
    
    options.progTol = .005;
           
    CPVector = reshape(obj.CP', numel(obj.CP), 1);
     
    CPVector = L1General2_BBSG(@obj.FirstOrderSolverInterface, CPVector, lambdaL1Deform*ones(size(CPVector)), options);
     
    obj.CP = reshape( CPVector, numel(CPVector)/obj.numSamples, obj.numSamples)'; 
   
%     % Deform and compute the criterion:
%     indRT = repmat(1:sizeDyadicRT, obj.numSamples, 1); % Interpolating sites
%     
%     DeformField = (obj.BsplDictDeform * obj.CP')';
%     
%     indKnot = indRT - DeformField; % Knots
%     
%     d = 3; % Order of the B-splines
%      
%     XCMS_Deform = zeros(sizeDyadicRT, obj.sizeMZ, obj.numSamples);
%     
%     % Initial deformation      
%     for j = 1 : obj.numSamples
%     
%         [IntBspl Bspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d,0);
% 
%         XCMS_Deform(:,:,j) = IntBspl * SeaMassCoef(:,:,j);
%          
%     end
%     
%     ctnOld = obj.ctnVarLogAnscombe2D(obj.NMASK .* XCMS_Deform, obj.RTWinCTN_LevelN);
%     
%     % A new normaliser: if improves the criterion more than the threshold
%     % then start the deformation again, otherwise go to the next level: 
%     NewNMASK = obj.NormaliseLS(XCMS_Deform, sizeNorm*2^(i-ResBegin+1)); 
%     
%     ctnNew = obj.ctnVarLogAnscombe2D(NewNMASK .* XCMS_Deform, obj.RTWinCTN_LevelN);
%         
%     if ctnOld - ctnNew < NThresh
%         
%         NFlag = 0;
%         
%     end
%     
     
end 
   
    
% Deformation Field at interpolation positions:
RT_Adjustment = (obj.BsplDictDeform * obj.CP')' ;
    
% Convert the index vector into minuts:
% RT_Adjustment = RT_Adjustment(:, obj.RTWinCTN_LevelN(1):obj.RTWinCTN_LevelN(2)) / (2^ResEnd);

end
