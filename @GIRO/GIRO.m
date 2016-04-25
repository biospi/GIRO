classdef GIRO < GIRO.MultiresRegister
   
    properties (Access = private)
                   
        lambdaL1Deform                         
                   
        CP = []      % Control Points for Deformation

        indCP        % Indicator of active CP 
        
        BsplDictDeform % Multiresolution B-spline deformation field
        
        lengthB_Deform
        
        lengthB_Norm
        
        sizeMZ
        
        OBJ_Normalise
        
        OBJ_CTN_LS_2D
        
        searchingStepSize
        
    end
    
    methods (Access = protected)
       
        % Generating deformed integrated B-spline dictionary for image representation:
        [IBspl, DIBspl] = BsplIntegral(OBJ_GIRO, t, x, d, NFactor)
        
        % Interface between register and solver defined at the Register class level:        
        [ctn, de_dc] = interf_register_solver(OBJ_GIRO, CP)
      
        % Generate a multi-resolution B-spline dictionary for deformation:
        [BsplDictDeform, normBuL2, numCoefs] = get_BsplDict_deform(OBJ_GIRO, numResLevels);

        OBJ_GIRO = searching_strategy(OBJ_GIRO, ~)
  
    end
    
    methods (Access = public)
        
        function OBJ_GIRO = GIRO(OBJ_Data, numLevels, lambdaL1Deform, lengthB_Deform, lengthB_Norm, solverDir, searchingStepSize)
           
            OBJ_GIRO@GIRO.MultiresRegister(OBJ_Data, numLevels);
            
            OBJ_GIRO.OBJ_Normalise = GIRO.Normalise(OBJ_Data);
            
            OBJ_GIRO.lambdaL1Deform = lambdaL1Deform;
            
            OBJ_GIRO.lengthB_Deform = lengthB_Deform;
            
            OBJ_GIRO.lengthB_Norm = lengthB_Norm;
            
            OBJ_GIRO.sizeMZ = size(OBJ_GIRO.Samples,2);
            
            OBJ_GIRO.OBJ_CTN_LS_2D = GIRO.CTN_LS_2D(OBJ_GIRO.Samples); 
            
            OBJ_GIRO.searchingStepSize = searchingStepSize;
            
            addpath(genpath(solverDir));
            
        end
        
        function Normaliser = get_Normaliser(OBJ_GIRO)
           
            Normaliser = OBJ_GIRO.OBJ_Normalise.get_Normaliser();
            
        end
             
    end
    
end
