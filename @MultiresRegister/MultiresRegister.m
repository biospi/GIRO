classdef MultiresRegister < GIRO.Register

    properties (Access = protected)
  
        Down_Samples % Downsampled images
        
        ResRT_Begin

        ResRT_End
        
        numLevels
        
        currentLevel
        
        sizeRT_currentLevel
        
        sizeDyadicRT_currentLevel
        
        indRT_Start_currentLevel
           
        LP_Taps = [.0417 .4583 .4583 .0417]
        
    end
    
    methods (Access = public)
        
        function OBJ_MultiresRegister = MultiresRegister(OBJ_Data, numLevels)
            
            OBJ_MultiresRegister@GIRO.Register(OBJ_Data);
            
            OBJ_MultiresRegister.numLevels = numLevels;
            
            OBJ_MultiresRegister.ResRT_End = OBJ_Data.get_levelsRT();
                        
            OBJ_MultiresRegister.ResRT_Begin = OBJ_MultiresRegister.ResRT_End - numLevels + 1;
            
            OBJ_MultiresRegister.Down_Samples = cell(OBJ_MultiresRegister.numLevels, 1);
            
            Samples = OBJ_MultiresRegister.Samples;
            
            for i = 1 : OBJ_MultiresRegister.numLevels
               
                % Lowpass filtering using the given LP_Taps:
                OBJ_MultiresRegister.Down_Samples(end - i + 1) = {Samples};
                
                Down_Samples = zeros( size(Samples(1:2:end,:,:),1), size(Samples,2), OBJ_MultiresRegister.numSamples );
                
                for j = 1 : OBJ_MultiresRegister.numSamples
                
                    Samples(:,:,j) = filter(OBJ_MultiresRegister.LP_Taps, 1, Samples(:,:,j), [], 1);
                
                    Down_Samples(:,:,j) = Samples(1:2:end,:,j);
                    
                end
                
                % Downsampling:
                Samples = Down_Samples;
                               
            end
      
            
        end
        
        
        function OBJ_MultiresRegister = start_multires_register(OBJ_MultiresRegister, SearchingParam)
            
            
            % Establish multi-resolution representation:
            for i = OBJ_MultiresRegister.ResRT_Begin : OBJ_MultiresRegister.ResRT_End
                
                OBJ_MultiresRegister.currentLevel = i;
                  
                % Call searching_strategy to do iterations within each resolution level:
                OBJ_MultiresRegister.searching_strategy(SearchingParam);
        
                
            end
            
        end
      
        
    end


end 







