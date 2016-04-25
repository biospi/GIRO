function OBJ_GIRO = searching_strategy(OBJ_GIRO, ~)                         
          
% Get the downsampled samples:
Down_Samples = cell2mat(OBJ_GIRO.Down_Samples(OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin + 1));

OBJ_GIRO.sizeRT_currentLevel = size(Down_Samples,1);

OBJ_GIRO.sizeDyadicRT_currentLevel = pow2( ceil(log2(OBJ_GIRO.sizeRT_currentLevel)));

OBJ_GIRO.indRT_Start_currentLevel = round( (OBJ_GIRO.sizeDyadicRT_currentLevel - OBJ_GIRO.sizeRT_currentLevel) / 2);

% Establish the level i multiresolution deformation dictionary:
[OBJ_GIRO.BsplDictDeform, normBL2, numCoefs] = OBJ_GIRO.get_BsplDict_deform(OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin + 1);

% Attach the new control points to the control points from previous levels:  
if (OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin) >= 1      

    indLevels = [0 cumsum(numCoefs(2:end))];   
    
    for j = 1 : (OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin) 
    
        OBJ_GIRO.CP(:, indLevels(j) + 1 : indLevels(j+1)) = OBJ_GIRO.CP(:, indLevels(j) + 1 : indLevels(j+1)) * 2 * normBL2(j+1) / normBL2(j);
        
    end
    
end

numCP_CurrentLevel = size(OBJ_GIRO.BsplDictDeform,2) - size(OBJ_GIRO.CP,2);

% Taken as initial of the next level:
OBJ_GIRO.CP = [zeros(OBJ_GIRO.numSamples, numCP_CurrentLevel), OBJ_GIRO.CP]; 

OBJ_GIRO.indCP = ones(size(OBJ_GIRO.CP));

% Initial Deformation:
indRT = repmat(1:OBJ_GIRO.sizeDyadicRT_currentLevel, OBJ_GIRO.numSamples, 1); % Interpolating sites

OBJ_GIRO.RT_Adjustment = (OBJ_GIRO.BsplDictDeform * OBJ_GIRO.CP')';

indKnot = indRT - OBJ_GIRO.RT_Adjustment; % Knots

d = 3; % Order of the B-splines

OBJ_GIRO.DeformedSamples = zeros(size(Down_Samples));

% Initial deformation      
for j = 1 : OBJ_GIRO.numSamples

    [IntBspl, DIBspl] = OBJ_GIRO.BsplIntegral(indKnot(j,:), indRT(j,:), d, 0);

    OBJ_GIRO.DeformedSamples(:,:,j) = IntBspl * Down_Samples(:,:,j);
    
end

% Normalisation:
OBJ_GIRO.OBJ_Normalise = OBJ_GIRO.OBJ_Normalise.update_normalise(OBJ_GIRO.DeformedSamples);
OBJ_GIRO.OBJ_Normalise = OBJ_GIRO.OBJ_Normalise.normalise_LS(OBJ_GIRO.lengthB_Norm*2^(OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin + 1)); %sizeNorm*2^(i+1)

% L1 regularised deform:
% Call the solver:
options.verbose = 0;

options.mode = 1;  

options.order = -1;

options.progTol = .1;

options.maxIter = 80;

% L1 regularised deformation:
CPVector = reshape(OBJ_GIRO.CP', numel(OBJ_GIRO.CP), 1);

% CPVector = L1General2_OWL(@OBJ_GIRO.interf_L1General_1stOrder_deform, CPVector, OBJ_GIRO.lambdaL1Deform * ones(size(CPVector)), options);
CPVector = L1General2_OWL(@OBJ_GIRO.interf_register_solver, CPVector, OBJ_GIRO.lambdaL1Deform * ones(size(CPVector)), options);
 
% LS deform:
OBJ_GIRO.CP = reshape( CPVector, numel(CPVector)/OBJ_GIRO.numSamples, OBJ_GIRO.numSamples)'; 

OBJ_GIRO.indCP = (OBJ_GIRO.CP ~= 0);

CPVector = L1General2_OWL(@OBJ_GIRO.interf_register_solver, CPVector, zeros(size(CPVector)), options);

OBJ_GIRO.CP = reshape( CPVector, numel(CPVector)/OBJ_GIRO.numSamples, OBJ_GIRO.numSamples)'; 

OBJ_GIRO.RT_Adjustment = -(OBJ_GIRO.BsplDictDeform * OBJ_GIRO.CP')' ;

OBJ_GIRO.RT_Adjustment = OBJ_GIRO.RT_Adjustment(:, OBJ_GIRO.indRT_Start_currentLevel : (OBJ_GIRO.sizeRT_currentLevel+OBJ_GIRO.indRT_Start_currentLevel-1)) * OBJ_GIRO.ResRT_Sec * 2^(OBJ_GIRO.ResRT_End - OBJ_GIRO.currentLevel);

OBJ_GIRO.OBJ_Normalise.normalise_median();

end
 