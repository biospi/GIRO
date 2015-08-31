function [ctn, gradCtn] = get_ctn_var_LogAnscombe_2D(OBJ_GIRO,Samples, RTWin) 
% Criterion for Registration: Variances
%
% function [ctn gradVar4Norm gradVar4Reg HVar] = ctnVarLogAnscombe2D(Samples) 
%
% This part of program receives the argument Samples and makes an
% evaluation of the minimisation criterion/gradient/hessian for the 
% registration.
%
% Samples now is a 3D matrix with the first two dim image and the third 
% dim sample numbers;
%
% RTWin is the retention time window which contains useful peptide signals
% 
Samples(Samples< 50) = 0; 

% Excluding zeros in all the samples:
zeroMask = ~(sum(Samples == 0, 3) ~= 0);

for i = 1 : OBJ_GIRO.numSamples
    
    tmp = Samples(:,:,i);
    
    tmp = tmp .* zeroMask;
    
    Samples(:,:,i) = tmp;
    
end

% Set the values outside retention window to be zero: excluded from
% criterion of deformation.
%
Samples([1:(RTWin(1)-1) (RTWin(end)+1):end], : , :) = 0;

Samples = log(Anscombe(Samples));

meanSamples = mean(Samples,3);

meanSamples = repmat(meanSamples,[1,1,OBJ_GIRO.numSamples]);

ctn = sum(sum(std( Samples , 0, 3))); % unconstrained criterion
 
if isnan(ctn)
    
    input('Criterion is NaN')
    
end

% ctn = sum(std( ParaSamples(1:numSamples,:) , 0, 1)) + LambdaSmooth * sum(Smooth(:));

gradCtn =  (Samples - meanSamples); % Other images will affect the Jacobian through the mean;
        
end
 