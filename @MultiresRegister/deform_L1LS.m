function [OBJ_GIRO, RT, RT_Adjustment, SamplesDeformed] = deform_L1LS(OBJ_GIRO,  lambdaL1Deform)

% This is an experimental version of the deformation algorithm, in the
% sense:
% 1. Data should be read and passed to the algorithm through another
% function, but in this experiment the data read algorithm is within the
% deformation algorithm. 
              
sizeNorm = 4; % Size between two knots of B-spline basis in LS normalisation 

BsplIntegMZ = [.0417 .4583 .4583 .0417];
                    
OBJ_GIRO.lambdaL1Deform = lambdaL1Deform;

% Retention window used for computing the criterion, meaningful retention
% window at each level is derived from this interval:

for i = OBJ_GIRO.SeaMassResRT_Begin : OBJ_GIRO.SeaMassResRT_End
    
    offsetMZ = zeros(1,OBJ_GIRO.numSamples);
    
    offsetRT = zeros(1,OBJ_GIRO.numSamples);
    
    cd([OBJ_GIRO.folderName num2str(i)])
    
    if i ~= 5
    
    % Load in the SeaMass coefs:    
    for j = 1 : OBJ_GIRO.numSamples
                    
        datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(1)) '/%d/%d/%d/%d/0/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
            
        offset = double(h5readatt(cell2mat(OBJ_GIRO.fileName(j)),datasetname, 'Offset'));
        
        offsetMZ(j) = offset(1);
        
        offsetRT(j) = offset(2)+1;
        
    end  

    else
        
    % Load in the offset of SeaMass coefs:
    for j = 1 : OBJ_GIRO.numSamples
                    
        datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(1)) '/%d/%d/%d/%d/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
       
        offset = double(h5readatt(cell2mat(OBJ_GIRO.fileName(j)),datasetname, 'Offset'));
        
        offsetMZ(j) = offset(1);
        
        offsetRT(j) = offset(2)+1;
        
    end  
    
    end

    OBJ_GIRO.offsetRT = max(offsetRT);

    offsetRT = OBJ_GIRO.offsetRT - offsetRT + 1;

    OBJ_GIRO.offsetMZ = max(offsetMZ); 

    offsetMZ = OBJ_GIRO.offsetMZ - offsetMZ + 1;

    RTWin_i = round(OBJ_GIRO.RTWin / 2^(OBJ_GIRO.SeaMassResRT_End - i)) + 1;
    
    sizeRT_i = RTWin_i(2) - RTWin_i(1) + 1;
    
    OBJ_GIRO.sizeDyadicRT = pow2(ceil(log2(sizeRT_i)));    
        
    OBJ_GIRO.RTWinCTN_LevelN(1) = round((OBJ_GIRO.sizeDyadicRT - sizeRT_i) / 2);
    
    OBJ_GIRO.RTWinCTN_LevelN(2) = OBJ_GIRO.RTWinCTN_LevelN(1) + sizeRT_i - 1;
        
    % Allocate a dyadically sized RT for the samples:
    SeaMassCoef = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples);
    
    indMZ = zeros(OBJ_GIRO.numChannels + 1, 1);
    
    for j = 1 : OBJ_GIRO.numSamples

        for k = 1 : OBJ_GIRO.numChannels
        
        cd([OBJ_GIRO.folderName num2str(i)])
        
        MZWin_k = cell2mat(OBJ_GIRO.MZWin(k));
        
        indMZ(k+1) = indMZ(k) + MZWin_k(2) - MZWin_k(1) + 1;
        
        if i ~= 5
        
        datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(k)) '/%d/%d/%d/%d/0/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);

        else
            
            datasetname = sprintf(['/' cell2mat(OBJ_GIRO.channelName(k)) '/%d/%d/%d/%d/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);

        end
        
        SeaMassCoef(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2), (indMZ(k)+1) : indMZ(k+1) , j) = double(h5read(cell2mat(OBJ_GIRO.fileName(j)),datasetname,[offsetMZ(j) offsetRT(j)],[MZWin_k(2)-MZWin_k(1) + 1 sizeRT_i])');
        
        end
        
        % MZ direction B-spline interpolation using a 4-tap FIR filter:
        SeaMassCoef(:,:,j) = filter(BsplIntegMZ, 1, SeaMassCoef(:,:,j), [], 2);
       
    end
        
    OBJ_GIRO.SeaMassCoef = SeaMassCoef;
    
    %% Establish the level i multiresolution deformation dictionary:
    [OBJ_GIRO.BsplDictDeform, normBuL2, numCoefs] = OBJ_GIRO.get_BsplDict_deform(i - OBJ_GIRO.SeaMassResRT_Begin + 1);
              
    % and the corresponding control points:
    % Attach the new control points to the control points from 
    % previous levels:  
    if (i - OBJ_GIRO.SeaMassResRT_Begin) >= 1      
    
        indLevels = [0 cumsum(numCoefs(2:end))];   
                 
    for j = 1 : (i-OBJ_GIRO.SeaMassResRT_Begin) 
        
        OBJ_GIRO.CP(:, indLevels(j) + 1 : indLevels(j+1)) = OBJ_GIRO.CP(:, indLevels(j) + 1 : indLevels(j+1)) * 2 * normBuL2(j+1) / normBuL2(j);
        
    end
    
    end 
     
    numCP_Level_i = size(OBJ_GIRO.BsplDictDeform,2) - size(OBJ_GIRO.CP,2);
    
    % Taken as initial of the next level:
    OBJ_GIRO.CP = [zeros(OBJ_GIRO.numSamples, numCP_Level_i), OBJ_GIRO.CP]; 

    OBJ_GIRO.indCP = ones(size(OBJ_GIRO.CP));
    
    %% Construct the integrated 3rd order B-spline basis matrix and its 
    %  derivative matrix: the 3rd order B-spline basis matrix
    indRT = repmat(1:OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.numSamples, 1); % Interpolating sites
    
    DeformField = (OBJ_GIRO.BsplDictDeform * OBJ_GIRO.CP')';
    
    indKnot = indRT - DeformField; % Knots
    
    d = 3; % Order of the B-splines
     
    SamplesDeformed = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples);
    
    % Initial deformation      
    for j = 1 : OBJ_GIRO.numSamples
    
        [IntBspl, DIBspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d,0);

        SamplesDeformed(:,:,j) = IntBspl * SeaMassCoef(:,:,j);
         
    end
    
     
   
    %% LS Normalising:
    OBJ_GIRO.NMASK = OBJ_GIRO.normalise_LS(SamplesDeformed, sizeNorm*2^(i-OBJ_GIRO.SeaMassResRT_Begin+1)); %sizeNorm*2^(i+1)
         
    %% Deforming:
        
%   % Call the solver:
    options.verbose = 0;
   
    options.mode = 1;  
  
    options.order = -1;
    
    options.progTol = .1;
           
    options.maxIter = 80;
    
    % L1 regularised deformation:
    CPVector = reshape(OBJ_GIRO.CP', numel(OBJ_GIRO.CP), 1);
     
    CPVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_deform, CPVector, OBJ_GIRO.lambdaL1Deform * ones(size(CPVector)), options);
        
    % Least square deformation using L1 selected control points:
    CP = reshape( CPVector, numel(CPVector)/OBJ_GIRO.numSamples, OBJ_GIRO.numSamples)'; 
    
    OBJ_GIRO.indCP = (CP ~= 0);
    
 %   CPVector = reshape(OBJ_GIRO.CP', numel(OBJ_GIRO.CP), 1);
    
    CPVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_deform, CPVector, zeros(size(CPVector)), options);
     
    OBJ_GIRO.CP = reshape( CPVector, numel(CPVector)/OBJ_GIRO.numSamples, OBJ_GIRO.numSamples)'; 
       
end 

% Deformation Field at interpolation positions:
RT_Adjustment = (OBJ_GIRO.BsplDictDeform * OBJ_GIRO.CP')' ;

indRT = repmat(1:OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.numSamples, 1); % Interpolating sites

indKnot = indRT - RT_Adjustment; % Knots

d = 3; % Order of the B-splines

SamplesDeformed = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples);

% Initial deformation
for j = 1 : OBJ_GIRO.numSamples

    [IntBspl, DIBspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d,0);

    SamplesDeformed(:,:,j) = IntBspl * SeaMassCoef(:,:,j);

end

RT = (OBJ_GIRO.offsetRT : (OBJ_GIRO.offsetRT + OBJ_GIRO.sizeRT - 1)) / (2^OBJ_GIRO.SeaMassResRT_End);

OBJ_GIRO.RT = RT;

% Convert the index vector into minuts:
RT_Adjustment = -RT_Adjustment(:, OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2)) / (2^OBJ_GIRO.SeaMassResRT_End);

SamplesDeformed = SamplesDeformed(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2),:,:);

OBJ_GIRO.SamplesDeformed = SamplesDeformed;

end 



