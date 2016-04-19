classdef (Abstract) CTN
    
    properties (Access =  private)
        
        Criterion
        
        gradCTN
        
    end
        
    methods (Access = public)
   
        function OBJ_CTN = OBJ_CTN(Samples, WinCTN)
                      
            [OBJ_CTN.Criterion, OBJ_CTN.gradCTN] = get_Criterion(Samples, WinCTN);
            
        end
        
        [Criterion, gradCTN] = get_Criterion(Samples, WinCTN);
        
    end
    
end
    