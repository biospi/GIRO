classdef CTN_LS_2D < CTN
    
    methods (Access = public)
        
        function [Criterion, gradCTN] = get_Criterion(Samples, WinCTN)
            
            zeroMask = ~(sum(Samples == 0, 3) ~= 0);
            
            for i = 1 : numSamples
                
                tmp = Samples(:,:,i);
                
                tmp = tmp .* zeroMask;
                
                Samples(:,:,i) = tmp;
                
            end
            
            
            % Set the values outside retention window to be zero: excluded
            % from criterion of deformation.

            Samples([1:(WinCTN(1)-1) (WinCTN(end)+1):end], : , :) = 0;

            Samples = log(Anscombe(Samples));

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