classdef Data
    
    properties (Access = protected)
        
        % LC-MS data related parameters:     
        folderName
        
        fileName
        
        Samples
        
        numSamples
        
        sizeRT
                                        
        sizeDyadicRT
        
        indRT_Start
        
        sizeMZ
        
        levelRT % log2(sizeRT)
        
        levelMZ % log2(sizeMZ)
        
        RT
        
        MZ

    end
     
    methods (Access = public)
                
        function OBJ_Data = Data(folderName, fileName)
            
            OBJ_Data.folderName = folderName;
        
            OBJ_Data.fileName = fileName;
                        
            OBJ_Data.numSamples = length(fileName);
            
        end
                 
        function Samples = get_Samples(OBJ_Data)
     
            Samples = OBJ_Data.Samples;
                        
        end
        
        function numSamples = get_numSamples(OBJ_Data)
            
            numSamples = OBJ_Data.numSamples;
            
        end
        
        function levelsRT = get_levelsRT(OBJ_Data)
            
            levelsRT = log2(OBJ_Data.sizeDyadicRT);
            
        end
        
        function ResRT_Sec = get_ResRT_Sec(OBJ_Data)
            
            ResRT_Sec = OBJ_Data.RT(2) - OBJ_Data.RT(1); % For uniformly sampled data
            
        end
          
        function [RT, sizeRT, sizeDyadicRT, indRT_Start] = get_RT_info(OBJ_Data)
        
            RT = OBJ_Data.RT;
            
            sizeRT = OBJ_Data.sizeRT;
            
            sizeDyadicRT = OBJ_Data.sizeDyadicRT;
            
            indRT_Start = OBJ_Data.indRT_Start;
            
        end
            
        function MZ = get_MZ(OBJ_Data)
        
            MZ = OBJ_Data.MZ;
            
        end
            
        function dispSample(OBJ_Data, indFig, indSample)
            
            figure(indFig); imagesc(OBJ_Data.MZ, OBJ_Data.RT, OBJ_Data(:,:,indSample))
            
        end
        
        dispOverlay(OBJ_Data)
        
        dispGroupOverlay(OBJ_Data)
               
        function [SampleLowResMZ, sizeMZ] = DyadicDownsampleMZ(OBJ_Data,SampleHighResMZ, numLevels, LP_Taps)
        
            for i = 1 : numLevels
               
                % Lowpass filtering using the given LP_Taps:
                SampleHighResMZ = filter(LP_Taps, 1, SampleHighResMZ, [], 2);
                
                % Downsampling:
                SampleHighResMZ = SampleHighResMZ(:, 1:2:end);
                               
            end
            
            SampleLowResMZ = SampleHighResMZ;
            
            sizeMZ = size(SampleLowResMZ, 2);
            
        end
        
                
        function SampleLowResRT = DyadicDownsampleRT(OBJ_Data, SampleHighResRT, numLevels, LP_Taps)
           
            for i = 1 : numLevels
               
                % Lowpass filtering using the given LP_Taps:
                SampleHighResRT = filter(LP_Taps, 1, SampleHighResRT, [], 1);
                
                % Downsampling:
                SampleHighResRT = SampleHighResRT(1:2:end, 1);
                               
            end
            
            SampleLowResRT = SampleHighResRT;
                        
        end
        
    end
               
end



