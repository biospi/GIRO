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
        
        numSamples = get_numSamples(OBJ_Data)
          
        [RT, sizeRT, sizeDyadicRT, indRT_Start] = get_RT_info(OBJ_Data)
        
        MZ = get_MZ(OBJ_Data)
        
        dispSample(OBJ_Data)
        
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
        
                
        function SampleLowResRT = DyadicDownsampleRT(OBJ_Data,SampleHighResRT, numLevels, LP_Taps)
           
            for i = 1 : numLevels
               
                % Lowpass filtering using the given LP_Taps:
                SampleHighResRT = filter(LP_Taps, 1, SampleHighResRT, [], 1);
                
                % Downsampling:
                SampleHighResRT = SampleHighResRT(1:2:end, 1);
                               
            end
            
            SampleLowResRT = SampleHighResMZ;
            
            
        end
        
    end
               
end



