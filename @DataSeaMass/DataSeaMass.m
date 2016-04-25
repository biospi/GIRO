classdef DataSeaMass < GIRO.Data
    
    properties (Access = private)
    
        numCommonChannels
        
        channelName
        
        SeaMassResMZ
        
        TargetResMZ
        
        SeaMassResRT
        
        SeaMassShrinkage
        
        SeaMassTolerance
                
        offsetRT_Sec
        
        offsetMZ_Thompson
                
        RTWin_Sec
        
        SecPerPixel
        
        MZWin_Thompson
        
        ThompsonPerPixel
                
        DownSampleFactorRT
        
        DownSampleFactorMZ
         
        % Internal constants                 
        BsplInteg = [.0417 .4583 .4583 .0417]
        
        C13_Constant = 1.00035
        
        Target_ResMZ_ThompsonPerPixel = 4
        
    end
        
    methods (Access = private)
        
        function OBJ_DataSeaMass = readData(OBJ_DataSeaMass)
                       
            offsetMZ = zeros(OBJ_DataSeaMass.numSamples, OBJ_DataSeaMass.numCommonChannels);
    
            offsetRT = zeros(OBJ_DataSeaMass.numSamples, OBJ_DataSeaMass.numCommonChannels);
    
            channelName = cell(OBJ_DataSeaMass.numSamples, OBJ_DataSeaMass.numCommonChannels);
            
            cd(OBJ_DataSeaMass.folderName)
            
            % Load in the offset of SeaMass coefs:
            for i = 1 : OBJ_DataSeaMass.numSamples
                               
                tmp = h5info(OBJ_DataSeaMass.fileName{i}, '/');
                tmp = tmp.Groups;
                numChannels = length(tmp);
                
                for j = 1 : OBJ_DataSeaMass.numCommonChannels
                    
                    for k = 1 : numChannels
                    
                    % Search the channel name for each of the common
                    % channels:
                    if ~isempty(regexp(tmp(k).Name, ['_' num2str(OBJ_DataSeaMass.channelName(j) )] ) ) 
                     
                        channelName(i,j) = {tmp(k).Name};
                        
                        datasetname = sprintf(['/' tmp(k).Name '/%d/%d/%d/%d/cs'], OBJ_DataSeaMass.SeaMassResMZ, OBJ_DataSeaMass.SeaMassResRT, OBJ_DataSeaMass.SeaMassShrinkage, OBJ_DataSeaMass.SeaMassTolerance);
       
                        offset = double(h5readatt(cell2mat(OBJ_DataSeaMass.fileName(i)), datasetname, 'Offset'));
       
                        offsetMZ(i, j) = offset(1);
       
                        offsetRT(i, j) = offset(2)+1;

                        break
                        
                    end
                    
                    end
                    
                end

            end

            % Compute the RT/MZ resolutions and coordinates for the data:            
            
            OBJ_DataSeaMass.offsetRT_Sec = max(offsetRT(:)) * OBJ_DataSeaMass.SecPerPixel;
            RTWin_Begin_Pixel = ceil( (OBJ_DataSeaMass.RTWin_Sec(1) - OBJ_DataSeaMass.offsetRT_Sec) / OBJ_DataSeaMass.SecPerPixel );
            offsetRT = (max(offsetRT(:)) - offsetRT + 1) + RTWin_Begin_Pixel;
            OBJ_DataSeaMass.RT = offsetRT(1)*OBJ_DataSeaMass.SecPerPixel : OBJ_DataSeaMass.SecPerPixel : (offsetRT(1) + (OBJ_DataSeaMass.sizeRT - 1) * OBJ_DataSeaMass.SecPerPixel);

            OBJ_DataSeaMass.sizeRT = size(OBJ_DataSeaMass.RT, 2);

            OBJ_DataSeaMass.offsetMZ_Thompson = max(offsetMZ(:)) * OBJ_DataSeaMass.ThompsonPerPixel; 
            MZWin_Begin_Pixel = ceil( (OBJ_DataSeaMass.MZWin_Thompson(1) - OBJ_DataSeaMass.offsetMZ_Thompson) / OBJ_DataSeaMass.ThompsonPerPixel );
            offsetMZ = max(offsetMZ(:)) - offsetMZ + MZWin_Begin_Pixel;
            OBJ_DataSeaMass.MZ = offsetMZ(1)*OBJ_DataSeaMass.ThompsonPerPixel : OBJ_DataSeaMass.ThompsonPerPixel : (offsetMZ(1) + (OBJ_DataSeaMass.sizeMZ - 1) * OBJ_DataSeaMass.ThompsonPerPixel);
          
            numMZ_Levels = ceil(log2(OBJ_DataSeaMass.Target_ResMZ_ThompsonPerPixel / OBJ_DataSeaMass.ThompsonPerPixel));
            
            sizeMZ = length(1 : 2^numMZ_Levels : OBJ_DataSeaMass.sizeMZ);
            
            % Read in the samples:
            OBJ_DataSeaMass.Samples = zeros(OBJ_DataSeaMass.sizeRT, sizeMZ * OBJ_DataSeaMass.numCommonChannels, OBJ_DataSeaMass.numSamples);
 
            for i = 1 : OBJ_DataSeaMass.numSamples
                
                for j = 1 : OBJ_DataSeaMass.numCommonChannels
                    
                    % Read in from smo
                    datasetname = sprintf(['/' channelName{i,j} '/%d/%d/%d/%d/cs'], OBJ_DataSeaMass.SeaMassResMZ, OBJ_DataSeaMass.SeaMassResRT, OBJ_DataSeaMass.SeaMassShrinkage, OBJ_DataSeaMass.SeaMassTolerance);

                    tmp = double(h5read(cell2mat(OBJ_DataSeaMass.fileName(i)), datasetname, [offsetMZ(i,j), offsetRT(i,j)], [OBJ_DataSeaMass.sizeMZ, OBJ_DataSeaMass.sizeRT]))';
                    
                    % mz downsampling:
                    if (numMZ_Levels > 0 )
                    
                    [tmp, sizeMZ] = OBJ_DataSeaMass.DyadicDownsampleMZ(tmp, numMZ_Levels, OBJ_DataSeaMass.BsplInteg);
                    
                    end
                    
                    % 
                    OBJ_DataSeaMass.Samples(:, (j-1)*sizeMZ+1 : j*sizeMZ, i) = tmp;                
                                        
                end
                
            end
            
            OBJ_DataSeaMass.sizeMZ = sizeMZ * OBJ_DataSeaMass.numCommonChannels;
            
            
        end
        
        function commonChannels = get_common_channels(OBJ_DataSeaMass)
                        
            
            channelName = cell(OBJ_DataSeaMass.numSamples,1);
            
            % Stage 1: get all the channel names in:
            for i = 1 : OBJ_DataSeaMass.numSamples
               
                tmp = h5info(OBJ_DataSeaMass.fileName{i}, '/');
                tmp = tmp.Groups;
                numChannels = length(tmp);
                
                precursorCenters = zeros(numChannels,1);
                
                for j = 1 : numChannels
                    
                    precursorCenterName = strsplit(tmp(j).Name, '_');
                    precursorCenters(j) = str2double(precursorCenterName(2));
                
                end
                
                channelName{i} = {precursorCenters};
                
            end
           
            % Stage 2: intersect all the channelName to get the common
            % channels:
            commonChannels = intersect(cell2mat(channelName{1}), cell2mat(channelName{2}));
            
            if OBJ_DataSeaMass.numSamples >= 3
                
                for i = 3 : OBJ_DataSeaMass.numSamples
                    
                    commonChannels = intersect(commonChannels, cell2mat(channelName{i}));
                    
                end
                
            end
               
            commonChannels = sort(commonChannels, 'ascend');
            
        end
            
    end
    
    methods (Access = public)
        
        function OBJ_DataSeaMass = DataSeaMass(folderName, fileName, SeaMassResRT, SeaMassResMZ, SeaMassShrinkage, SeaMassTolerance, RTWin_Sec, MZWin_Thompson)
        
            OBJ_DataSeaMass = OBJ_DataSeaMass@GIRO.Data(folderName, fileName);
                        
            OBJ_DataSeaMass.channelName = get_common_channels(OBJ_DataSeaMass);
            OBJ_DataSeaMass.numCommonChannels = length(OBJ_DataSeaMass.channelName);

            OBJ_DataSeaMass.SeaMassResRT = SeaMassResRT;
            OBJ_DataSeaMass.SeaMassResMZ = SeaMassResMZ;
            OBJ_DataSeaMass.SeaMassShrinkage = SeaMassShrinkage;
            OBJ_DataSeaMass.SeaMassTolerance = SeaMassTolerance;
             
            OBJ_DataSeaMass.RTWin_Sec = RTWin_Sec;      
            OBJ_DataSeaMass.SecPerPixel =  60 / (2^OBJ_DataSeaMass.SeaMassResRT);            
            OBJ_DataSeaMass.sizeRT = ceil( (RTWin_Sec(2) - RTWin_Sec(1)) / OBJ_DataSeaMass.SecPerPixel );
            OBJ_DataSeaMass.levelRT = ceil(log2(OBJ_DataSeaMass.sizeRT));
            OBJ_DataSeaMass.sizeDyadicRT = pow2(OBJ_DataSeaMass.levelRT); 
            OBJ_DataSeaMass.indRT_Start = ceil( (OBJ_DataSeaMass.sizeDyadicRT - OBJ_DataSeaMass.sizeRT) / 2);
            
            OBJ_DataSeaMass.MZWin_Thompson = MZWin_Thompson;
            OBJ_DataSeaMass.ThompsonPerPixel = OBJ_DataSeaMass.C13_Constant / (60 * (2^OBJ_DataSeaMass.SeaMassResMZ));            
            OBJ_DataSeaMass.sizeMZ = ceil( (MZWin_Thompson(2) - MZWin_Thompson(1)) / OBJ_DataSeaMass.ThompsonPerPixel );
            OBJ_DataSeaMass.levelMZ = floor(log2(OBJ_DataSeaMass.sizeMZ));
            
            OBJ_DataSeaMass = readData(OBJ_DataSeaMass);
            
        end
        

        
        
        
    end
    
end