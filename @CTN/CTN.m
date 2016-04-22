classdef CTN
    
    properties (Access =  private)
        
        Criterion
        
        gradCTN
        
    end
        
    methods (Access = public)
   
        function OBJ_CTN = CTN(Samples)
                      
            [OBJ_CTN.Criterion, OBJ_CTN.gradCTN] = OBJ_CTN.get_Criterion(Samples);
            
        end
        
        [Criterion, gradCTN] = get_Criterion(OBJ_CTN, Samples);
        
    end
    
end
    