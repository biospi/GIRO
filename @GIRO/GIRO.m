classdef GIRO
    
    properties (Access = private)
        
        % Data file related parameters:
        fileName
        
        folderName
        
        channelName
        
        SeaMassResMZ
                
        SeaMassResRT_Begin
        
        SeaMassResRT_End
                
        SeaMassShrinkage
        
        SeaMassTolerance 
        
        % LC-MS data related parameters:
        SeaMassCoef
        
        sizeRT
        
        sizeDyadicRT
         
        sizeMZ
        
        offsetRT
        
        offsetMZ
        
        NMASK
        
        RT
        
        MZ
        
        resoluteRT
        
        resoluteMZ
        
        RTWin % Retention window for loading in images
        
        MZWin % M/Z window for loading in images
        
        RTWinCTN_LevelN
        
        numSamples
        
        numChannels
       
        indChannel

        % Details of image representation and normalisation parameters
        % should be defined in different methods. 
                
        lambdaL1Deform
                             
        CP           % Control Points for Deformation

        indCP        % Indicator of CP selected by L1 
        
        BsplDictDeform % Multiresolution B-spline deformation field
        
        % Parameters for Normalisation:
       
        SamplesDeformed % Deformed images
        
        SamplesNormalised % Deformed, normalised images
       
        BsplDictNorm % Multi-resolution B-spline normalisation basis
       
        lambdaL1Norm
       
        CN % Coefficients for normalisation to be multiplied with BsplDictNorm
       
        indCN % Indicator of CN selected by L1 
        
        RefNorm % Target for normalisation
        
        TarNorm % Source for normalisation
       
    end
    
    methods (Access = private)
        
        % Multi-resolution image representation schemes:
        % Re-binning
        
        % Normalisation schemes:
        NMASK = normalise_LS(OBJ_GIRO, Samples, LengthBu);
                           
        [ctn, gradCtn] = get_ctn_var_LogAnscombe_2D(OBJ_GIRO,Samples, RTWin); 
           
        % Generate a multi-resolution B-spline dictionary for deformation:
        [BsplDictDeform, normBuL2, numCoefs] = get_BsplDict_deform(OBJ_GIRO, numResLevels);

        % Generate a multi-resolution B-spline dictionary for normalisation:
        [numBspl_1D, numBspl_2D, B] = get_BsplDict_norm(OBJ_GIRO, minScale1D, maxScale1D, minScale2D, maxScale2D, sizeMZ_Channel, normalisation);

    end
     
    methods (Access = public)
                
        function OBJ_GIRO = GIRO(folderName, fileName, channelName, SeaMassResMZ, SeaMassResRT_Begin, SeaMassResRT_End, SeaMassShrinkage, SeaMassTolerance, RTWin, MZWin)
            
            OBJ_GIRO.folderName = folderName;
        
            OBJ_GIRO.fileName = fileName;
            
            OBJ_GIRO.channelName = channelName;
            
            OBJ_GIRO.SeaMassResMZ = SeaMassResMZ;
            
            OBJ_GIRO.SeaMassResRT_Begin = SeaMassResRT_Begin;
            
            OBJ_GIRO.SeaMassResRT_End = SeaMassResRT_End;
            
            OBJ_GIRO.SeaMassShrinkage = SeaMassShrinkage;
            
            OBJ_GIRO.SeaMassTolerance = SeaMassTolerance;
            
            OBJ_GIRO.RTWin = RTWin;
            
            OBJ_GIRO.MZWin = MZWin;
            
            OBJ_GIRO.numSamples = length(fileName);
            
            OBJ_GIRO.numChannels = length(channelName);

            OBJ_GIRO.indChannel = ones(OBJ_GIRO.numChannels+1,1);
            
            OBJ_GIRO.sizeRT = RTWin(2)-RTWin(1)+1;
 
            OBJ_GIRO.sizeMZ = 0;
            
            OBJ_GIRO.CP = [];
            
            OBJ_GIRO.indCP = 1;
            
            for k = 1 : OBJ_GIRO.numChannels
            
                MZWin_k = cell2mat(MZWin(k));
                
                OBJ_GIRO.sizeMZ = OBJ_GIRO.sizeMZ + MZWin_k(2)-MZWin_k(1)+1;
                                
                OBJ_GIRO.indChannel(k+1) = OBJ_GIRO.sizeMZ + 1; 
                
            end
                        
         end
          
        % L1 constrained registration using the mass-rebinned samples from 
        % SeaMass. The standard general L1 solver is used.
        
        [OBJ_GIRO, RT, RT_Adjustment, Samples_Deformed] = deform_L1LS(OBJ_GIRO, lambdaL1Deform);
                
        [ctn, de_dc] = interf_L1General_1stOrder_deform(OBJ_GIRO, CP);
                       
        % L1 constrained normalisation using the mass-rebinned samples from
        % SeaMass.  The standard general L1 solver is used.
        
        [MZ_Adjustment, Samples_Normalised] = normalise_L1LS(OBJ_GIRO, lambdaL1Norm);
        
        [ctn, de_dc] = interf_L1General_1stOrder_norm(OBJ_GIRO, CN);
        
        
    end
               
end



