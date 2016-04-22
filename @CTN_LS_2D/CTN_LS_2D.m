classdef CTN_LS_2D < GIRO.CTN
    
    methods (Access = public)
        
        function OBJ_CTN_LS_2D = CTN_LS_2D(Samples)
           
            OBJ_CTN_LS_2D@GIRO.CTN(Samples)
            
        end
        
        function [Criterion, gradCTN] = get_Criterion(OBJ_CTN_LS_2D, Samples)
            
            zeroMask = ~(sum(Samples == 0, 3) ~= 0);

            numSamples = size(Samples,3);
            
            for i = 1 : numSamples
                
                tmp = Samples(:,:,i);
                
                tmp = tmp .* zeroMask;
                
                Samples(:,:,i) = tmp;
                
            end
            
            
            % Set the values outside retention window to be zero: excluded
            % from criterion of deformation.
          
            Samples = log(GIRO.Anscombe(Samples));

            meanSamples = mean(Samples,3);

            meanSamples = repmat(meanSamples,[1,1,numSamples]);

            % unconstrained criterion
            Criterion = sum(sum(std( Samples , 0, 3))); 
 
            if isnan(Criterion)
    
                input('Criterion is NaN')
    
            end

            gradCTN =  (Samples - meanSamples); % Other images will affect the Jacobian through the mean;
            
        end 
        
    end
    
end