% function [XCMS_GIRO RT_Adjustment] = L1_LS_BsplPyrDeform(obj)

% This is an experimental version of the deformation algorithm, in the
% sense:
% 1. Data should be read and passed to the algorithm through another
% function, but in this experiment the data read algorithm is within the
% deformation algorithm. 
addpath D:\Data\David_Julian\Orbitrap1\SeaMass2D\MultiResolutions

fileName = ['20120531_DN1a_1_102p4_1_01_1.h5';.../
            '20120531_DN8a_1_102p4_1_01_1.h5'];

estimatedMaxDeform = 55;

numResLevels = floor(log2(estimatedMaxDeform));
        
numSamples = 2;        
        
numResLevels = length(fileName);        

sizeNorm = 4; % Size between two knots of B-spline basis in normalisation

BsplIntegMZ = [.0417 .4583 .4583 .0417];

% Retention window to load in from the .h5 file:
RTWin = [1 6140];

% Retention window used for computing the criterion, meaningful retention
% window at each level is derived from this interval:
RTWinCTN = [1035 5600];

% MZ window to load in from the .h5 file:
MZWin = [17 624];

sizeRT = 2^ceil(log(RTWin(end)-RTWin(1)));

sizeMZ = MZWin(2)-MZWin(1)+1;

% Allocate a dyadically sized RT size for the samples:
SeaMassCoef = zeros(sizeRT, sizeMZ, numSamples);

for j = 1 : numSamples
    
    % Read in the SeaMass coefs:
    SeaMassCoef(:,:,i) = double(h5read(fileName_i(j,:),'/Image',[RTWin(1) MZWin(1)], [RTWin(2)-RTWin(1) MZWin(2)-MZWin(1)] ) );
    
    % MZ direction B-spline interpolation using a 4-tap FIR filter:
    SeaMassCoef(:,:,j) = filter(BsplIntegMZ, 1, SeaMassCoef(:,:,j), [], 2);
    
end
    
SeaMassCoef = SeaMassCoef(:, (1+length(BsplIntegMZ)):(end-length(BsplIntegMZ)), :); % Remove the boundary

% Adjust the MZ parameters accordingly:
MZWin = [MZWin(1)+length(BsplIntegMZ) MZWin(2)-length(BsplIntegMZ)];
sizeMZ = MZWin(2)-MZWin(1)+1;

%% Iterative normalisation-deformation: 
NFlag = 1; % Normalising flag: 1 if 

NRThresh = 0;
 
lambdaL1Deform = .01;  

% In this version of GIRO, the B-spline deformation field can be defined
% upfrount. Only the coefs passing the threshold will be updated.
%[obj.BsplDictDeform BsplDictDeformInterp normBuL2 numCoefs] = BsplDeformDict(obj, obj.sizeRT, resoluteInterp, obj.numResLevels);

%obj.CP = zeros(obj.numSamples, size(obj.BsplDictDeform,1)); 
 
%%

for i = 1 : numResLevels
   
    % Establish the level i integrated B-spline basis for image
    % representation:
    
    
    % Establish the level i multiresolution deformation dictionary:
        
    % and the corresponding control points:
    
        
    %% Normalising:
    NMASK = NormaliseLS(XCMS, sizeNorm*2^(i+1));
   
    XCMS_Normalised = XCMS_Filtered .* NMASK; % This normalised XCMS will be used to run the deformation until no further updates can be done.   
    
    %% Deforming:
       
    % Evaluate the criterion and the derivative of image intensity w.r.t.
    % the criterion: de_dI
    RTWin_i = [floor(RTWin(1)/(2^(i-1))) ceil(RTWin(2)/(2^(i-1)))];
    
    [ctn de_dI] = ctnVarLogAnscombe2D(XCMS_Normalised, RTWin_i);
    
    % Establish the level i integrated B-spline basis for image
    % representation:
    
    % Establish the level i multiresolution deformation dictionary:
    [BsplDictDeform normBuL2 numCoefs] = BsplDeformDict(sizeRT_i, numCurrentResLevels);
          
    cumNormBuL2(numResLevels - i + 1) = {normBuL2};
 
    % and the corresponding control points:
       % Attach the new control points to the control points from 
    % previous levels: 
    if numCurrentResLevels >= 2  
        
        indLevels = [0 cumsum(numCoefs(2:end))];   
        
        % Get the L2 norm of previous level, 
        PreviousL2Norm = cell2mat(cumNormBuL2(numCurrentResLevels-1));
         
    for j = 1 : (numCurrentResLevels-1)
                 
        obj.CP(:, indLevels(j) + 1 : indLevels(j+1)) = obj.CP(:, indLevels(j) + 1 : indLevels(j+1)) * 2 * normBuL2(j+1) / PreviousL2Norm(j);
        
    end
    
    end 
    
    % Taken as initial of the next level:
    CP = [zeros(numSamples, size(BsplDictDeform,2) - size(CP,2)), obj.CP]; 
 
    % Call the solver:
    options.verbose = 0;
   
    options.mode = 1; 
  
    options.order = 1;
    
    options.progTol = .001;
           
    CPVector = reshape(obj.CP', numel(obj.CP), 1);
     
    CPVector = L1General2_BBSG(@obj.FirstOrderSolverInterface, CPVector, lambdaL1Deform*ones(size(CPVector)), options);
     
    obj.CP = reshape( CPVector, numel(obj.CP)/obj.numSamples, obj.numSamples)'; 
         
end
   
% 
% 
% Deformation Field at interpolation positions:
% RT_Adjustment = (BsplDictDeformInterp * CP')' ;
% 
% indOrg = repmat(.5/resoluteInterp : 1/resoluteInterp : (sizeRT- .5/resoluteInterp),obj.numSamples,1) - RT_Adjustment;
% 
% indNew = 1 : sizeRT;
% 
% XCMS_GIRO = obj.Interp1D_MassDist(obj.XCMS_RAW, indOrg, indNew, resoluteInterp);
