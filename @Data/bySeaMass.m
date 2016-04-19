function [RT, MZ, Intensity] = bySeaMass(channelName, fileName, SeaMassResRT, SeaMassResMZ, SeaMassShrinkage, SeaMassTolerance)    


BsplIntegMZ = [.0417 .4583 .4583 .0417];

offsetMZ = zeros(1,OBJ_GIRO.numSamples);
    
    offsetRT = zeros(1,OBJ_GIRO.numSamples);
    
    cd([OBJ_GIRO.folderName num2str(i)])
    
    if i >= 0
    
    % Load in the SeaMass coefs:    
    for j = 1 : OBJ_GIRO.numSamples
                    
        datasetname = sprintf(['/' cell2mat(channelName(1)) '/%d/%d/%d/%d/0/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
%        datasetname = sprintf(['/' cell2mat(channelName(1)) '/%d/%d/%d/%d/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
            
        offset = double(h5readatt(cell2mat(fileName(j)),datasetname, 'Offset'));
        
        offsetMZ(j) = offset(1);
        
        offsetRT(j) = offset(2)+1;
        
    end  

    else
        
    % Load in the offset of SeaMass coefs:
   for j = 1 : OBJ_GIRO.numSamples
                    
       datasetname = sprintf(['/' cell2mat(channelName(1)) '/%d/%d/%d/%d/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);
       
       offset = double(h5readatt(cell2mat(fileName(j)),datasetname, 'Offset'));
        
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
    
    OBJ_GIRO.sizeDyadicRT = 2*pow2(ceil(log2(sizeRT_i)));    
        
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
        
%        if i ~= 5
        
        datasetname = sprintf(['/' cell2mat(channelName(k)) '/%d/%d/%d/%d/0/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);

%        else
            
%            datasetname = sprintf(['/' cell2mat(channelName(k)) '/%d/%d/%d/%d/cs'], OBJ_GIRO.SeaMassResMZ, i, OBJ_GIRO.SeaMassShrinkage, OBJ_GIRO.SeaMassTolerance);

%        end
        
        SeaMassCoef(OBJ_GIRO.RTWinCTN_LevelN(1):OBJ_GIRO.RTWinCTN_LevelN(2), (indMZ(k)+1) : indMZ(k+1) , j) = double(h5read(cell2mat(fileName(j)),datasetname,[offsetMZ(j)+MZWin_k(1)-1 offsetRT(j)+RTWin_i(1)-1],[MZWin_k(2)-MZWin_k(1) + 1 sizeRT_i])');
        
        end
        
        % MZ direction B-spline interpolation using a 4-tap FIR filter:
        SeaMassCoef(:,:,j) = filter(BsplIntegMZ, 1, SeaMassCoef(:,:,j), [], 2);
       
    end
    
end
