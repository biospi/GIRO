function [XCMS_REP dI_dg] = LS_XCMS_Representation(obj, Samples, LengthBu)
% Least square image representation:            
 
sizeRT = size(Samples,1);

sizeMZ = size(Samples,2);

% B-spline basis table for image interpolation: a much finer resolution
%LengthBu = 4;

numCoefs = sizeRT/LengthBu;

% Allocate memory for sparse deforma basis:
BsplDictRep = spalloc(sizeRT, numCoefs, numCoefs*4*LengthBu);

BsplGradDict= spalloc(sizeRT, numCoefs, numCoefs*4*LengthBu);

% LengthBuInterp = LengthBu*resoluteInterp;

% Allocate memory for dense B-spline coefs: 
BsplInterpCoef = zeros(numCoefs, sizeMZ, obj.numSamples);

XCMS_REP = zeros(sizeRT,sizeMZ,obj.numSamples);

dI_dg = zeros(sizeRT, sizeMZ, obj.numSamples);


% Dictionary for representation:
%u = (.5 / LengthBu ) : 1/LengthBu : (1- .5/LengthBu);
u = (0 : LengthBu)/LengthBu;
u = u(1:(end-1));

Bu3 = u.^3/6;
Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;
Bu0 = (1-u).^3/6;

Bu = [Bu3 Bu2 Bu1 Bu0];
 
%BsplDictRep( 1 : LengthBu,1) = Bu0 / norm(Bu0,1);  

BsplDictRep( 1 : 2*LengthBu, 1) = [Bu1 Bu0] / norm([Bu1 Bu0],1);

BsplDictRep( 1 : 3*LengthBu, 2) = [Bu2 Bu1 Bu0] / norm([Bu2 Bu1 Bu0],1);

%BsplDictRep( end - LengthBu + 1 : end,end) = Bu3 / norm(Bu3, 1);

%BsplDictRep( end - 2*LengthBu + 1 : end,end-1) = [Bu3 Bu2] / norm([Bu3 Bu2],1);

BsplDictRep( end - 3*LengthBu + 1 : end, end) = [Bu3 Bu2 Bu1] / norm([Bu3 Bu2 Bu1],1);

normBu1 = norm(Bu,1);

Bu = Bu/normBu1;

for k = 3 : numCoefs-1
        
    BsplDictRep((k-3) * LengthBu + 1 : (k+1) * LengthBu, k) = Bu;

end

% Gradient of B-splines:
dBu3 = u.^2/2;
dBu2 = -1.5*u.^2 + u + .5;
dBu1 = 1.5*u.^2 - 2*u;
dBu0 = -(1-u).^2/2;

dBu = [dBu3 dBu2 dBu1 dBu0] / normBu1;

%BsplGradDict( 1 : LengthBu,1) = dBu0 / normBu1;  

BsplGradDict( 1 : 2*LengthBu, 1) = [dBu1 dBu0] / normBu1;

BsplGradDict( 1 : 3*LengthBu, 2) = [dBu2 dBu1 dBu0] / normBu1;

%BsplGradDict( end - LengthBu + 1 : end,end) = dBu3 / normBu1;

%BsplGradDict( end - 2*LengthBu + 1 : end,end-1) = [dBu3 dBu2] / normBu1;

BsplGradDict( end - 3*LengthBu + 1 : end, end) = [dBu3 dBu2 dBu1] / normBu1;

for k = 3 : numCoefs-1
        
    BsplGradDict((k-3) * LengthBu + 1 : (k+1)*LengthBu,k) = dBu;

end


% Least square image representation:
for j = 1 : obj.numSamples
    
    % Unconstrained least square representation: may introduce negative coefs
    BsplInterpCoef(:,:,j) = BsplDictRep \ Samples(:,:,j);    

    % Non-negative least squares:
%     for k = 1 : sizeMZ
%         
%         BsplInterpCoef(:,k,j) = lsqnonneg(BsplDictRep, Samples(:,k,j));
%          
%     end
    
    XCMS_REP(:,:,j) = BsplDictRep * BsplInterpCoef(:,:,j);

    % the image gradient for dI_dg:
    dI_dg(:,:,j) = BsplGradDict * BsplInterpCoef(:,:,j);
 
end
    


