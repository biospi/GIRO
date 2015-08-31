function obj = LS_RigidDeform(obj, resL)

% The purpose of this function is to groupwise rigid register the images so
% that the variable offsetRT can be set before the B-spline deformation.

% Read in H5 files:
mz_resolution = -6;

shrinkage = -4;

tolerance = -9;

% SeaMassCoef = obj.SeaMassCoef.* obj.NMASK;
SeaMassCoef = obj.SeaMassCoef;

offsetRT = zeros(1,obj.numSamples);

offsetRT_Tmp = zeros(1,obj.numSamples);

deltaRT = 10; % Move up to 2 pixel for each iteration

deltaCtn = 10;

numIter = 1;

while (deltaCtn > 0) && (numIter < 120)
        
    % Find the derivative of the image deformation and integrate along MZ:
    dI_dg = zeros(obj.sizeDyadicRT, obj.sizeMZ, obj.numSamples);

    IntBspl = [0.0417 0.4583 0.4583 0.0417];

    DBspl = [0.1667 0.5000 -0.5000 -0.1667];
    
    % Deformation using B-splines:      
    for j = 1 : obj.numSamples

        SeaMassCoef(:,:,j) = filter(IntBspl,1,SeaMassCoef(:,:,j),[],1);
    
        dI_dg(:,:,j) = filter(DBspl,1, SeaMassCoef(:,:,j),[],1);
    
    end
  
    % Deforming according to CP:
      
    % B-spline image representation:
           
    % Criterion re-evaluation: 
    [ctn de_dI] = obj.ctnVarLogAnscombe2D(SeaMassCoef, obj.RTWinCTN_LevelN);  

    parfor j = 1 : obj.numSamples 
     
        % sub-gradient including the threshold:
        offsetRT_Tmp(j) = sum(sum(de_dI(:,:,j) .* dI_dg(:,:,j) ./ (obj.sizeMZ*(SeaMassCoef(:,:,j)+.375))));
   
    end 
    
    offsetRT = offsetRT + round(deltaRT * offsetRT_Tmp / max(offsetRT_Tmp))
  
    for j = 1 : obj.numSamples    
        % Check the new criterion:
        datasetname = sprintf('/1/%d/%d/%d/%d/0/cs', mz_resolution, resL, shrinkage, tolerance);
        
        SeaMassCoef(obj.RTWinCTN_LevelN(1):obj.RTWinCTN_LevelN(2),:,j) = double(h5read(obj.fileName(j,:),datasetname,[1+obj.offsetMZ(j) 1+offsetRT(j)],[obj.sizeMZ obj.sizeRT_i])');

        % MZ direction B-spline interpolation using a 4-tap FIR filter:
        SeaMassCoef(:,:,j) = filter(IntBspl,1,SeaMassCoef(:,:,j),[],2);
        
        SeaMassCoef(:,:,j) = filter(IntBspl,1,SeaMassCoef(:,:,j),[],1);
        
    end
 
%    SeaMassCoef = SeaMassCoef .* obj.NMASK;
    
    [ctnNew de_dI] = obj.ctnVarLogAnscombe2D(SeaMassCoef, obj.RTWinCTN_LevelN); 
    
    deltaCtn = ctn - ctnNew
    
    if deltaCtn > 0
        
        ctn =  ctnNew;
        
    end
    
end
    
    
    
    