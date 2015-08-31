classdef BsplNormDict 
   
    properties ( Access = private )
        
        VEC
        
    end
               
    methods ( Access = public )
%        
%         function objBsplNormDict = BsplNormDict(a,b)
%             
%             objBsplNormDict.MAT = a;
%             
%             objBsplNormDict.VEC = b;
%                                     
%         end
        
        function innerProd = MAT(objBsplNormDict, b)
            
            
            
            innerProd = magic(5) * b;
            
        end
                   
        function c = mtimes(a, b)
            
            c = a.MAT(b);
            
        end
        
    end
        
    
end