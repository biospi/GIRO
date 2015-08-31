function [ctn gradCtn] = ctnVar2D(obj,Samples) 
% Criterion for Registration: Variances
%
% function [ctn gradVar4Norm gradVar4Reg HVar] = ctnVar(Samples) 
%
% This part of program receives the argument Samples and makes an
% evaluation of the minimisation criterion/gradient/hessian for the 
% registration.
% Samples now is a 3D matrix with the first two dim image and the third 
% dim sample numbers;
meanSamples = mean(Samples,3);

meanSamples = repmat(meanSamples,[1,1,obj.numSamples]);

ctn = sum(sum(std( Samples , 0, 3))); % unconstrained criterion

if isnan(ctn)
    
    input('Criterion is NaN')
    
end

% ctn = sum(std( ParaSamples(1:numSamples,:) , 0, 1)) + LambdaSmooth * sum(Smooth(:));

gradCtn = 2 * (Samples - meanSamples); % Other images will affect the Jacobian through the mean;
        
end
 