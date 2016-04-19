classdef MultiresRegister < Register

    properties (Access = protected)
  
        Down_Samples % Downsampled images
        
        ResRT_Begin

        ResRT_End
        
        numLevels
        
        currentLevel
        
        RTWinCTN_LevelN    
        
    end
    
    methods (Access = public)
        
        function OBJ_MultiresRegister = MultiresRegister(OBJ_Data, numLevels)
            
            OBJ_MultiresRegister@Register(OBJ_Data);
            
            OBJ_MultiresRegister.numLevels = numLevels;
                        
        end
        
        
        function OBJ_MultiresRegister = start_multires_register(OBJ_MultiresRegister)
            
            % Establish multi-resolution representation:
            for i = OBJ_MultiresRegister.ResRT_Begin : OBJ_MultiresRegister.ResRT_End
            
                OBJ_MultiresRegister.currentLevel = i;
                
                % Get the multiresolution representation in place:
                
                
                % Call searching_strategy to do iterations within each resolution level:
                OBJ_MulitiresRegister.search_strategy();
        
%      OBJ_GIRO.channelName = channelName;
%             
%             OBJ_GIRO.SeaMassResMZ = SeaMassResMZ;
%             
%             OBJ_GIRO.SeaMassResRT_Begin = SeaMassResRT_Begin;
%             
%             OBJ_GIRO.SeaMassResRT_End = SeaMassResRT_End;
%             
%             OBJ_GIRO.SeaMassShrinkage = SeaMassShrinkage;
%             
%             OBJ_GIRO.SeaMassTolerance = SeaMassTolerance;
%             
%             OBJ_GIRO.RTWin = RTWin;
%             
%             OBJ_GIRO.MZWin = MZWin;
%             
%             OBJ_GIRO.numSamples = length(fileName);
%             
%             OBJ_GIRO.numChannels = length(channelName);
% 
%             OBJ_GIRO.indChannel = ones(OBJ_GIRO.numChannels+1,1);
%             
%             OBJ_GIRO.sizeRT = RTWin(2)-RTWin(1)+1;
%  
%             OBJ_GIRO.sizeMZ = 0;
%             
%             OBJ_GIRO.CP = [];
%             
%             OBJ_GIRO.indCP = 1;
%             
%             for k = 1 : OBJ_GIRO.numChannels
%             
%                 MZWin_k = cell2mat(MZWin(k));
%                 
%                 OBJ_GIRO.sizeMZ = OBJ_GIRO.sizeMZ + MZWin_k(2)-MZWin_k(1)+1;
%                                 
%                 OBJ_GIRO.indChannel(k+1) = OBJ_GIRO.sizeMZ + 1; 
%                 
%          
%                 
%                 
%                 
%                             
%             RTWin_i = round(OBJ_DataSeaMass.RTWin / 2^(OBJ_DataSeaMass.SeaMassResRT_End - i)) + 1;
%   
%             sizeRT_i = RTWin_i(2) - RTWin_i(1) + 1;
%   
%             OBJ_DataSeaMass.sizeDyadicRT = pow2(ceil(log2(sizeRT_i)));    
%   
%             OBJ_DataSeaMass.RTWinCTN_LevelN(1) = round((OBJ_DataSeaMass.sizeDyadicRT - sizeRT_i) / 2);
%   
%             OBJ_DataSeaMass.RTWinCTN_LevelN(2) = OBJ_DataSeaMass.RTWinCTN_LevelN(1) + sizeRT_i - 1;
%             
%             
%                 
%         % L1 constrained registration using the mass-rebinned samples from 
%         % SeaMass. The standard general L1 solver is used.
% %         
% %         [OBJ_GIRO, RT, RT_Adjustment, Samples_Deformed] = deform_L1LS(OBJ_GIRO, lambdaL1Deform);
% %                 
% %         [ctn, de_dc] = interf_L1General_1stOrder_deform(OBJ_GIRO, CP);
% %                        
% %         % L1 constrained normalisation using the mass-rebinned samples from
% %         % SeaMass.  The standard general L1 solver is used.
% %         
% %         [MZ_Adjustment, Samples_Normalised] = normalise_L1LS(OBJ_GIRO, lambdaL1Norm);
% %         
% %         [ctn, de_dc] = interf_L1General_1stOrder_norm(OBJ_GIRO, CN);
% 
% 
            end
            
        end
            
    end


end 







