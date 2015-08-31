function [BsplDict BsplInterpCoef BsplDictLocation] = BsplDict_1D(obj)

% This function is a private member method for the class LCMS: it
% initialise the dictionaries used for the LS image representations:
% 
% [BsplDictInv BsplDictLocation] = BsplDict_1D(obj)
%
% the output: BsplDictInv is the pseudo-inverse of the B-spline dictionary
% BsplGradDict is the gradient of the dictionary
% BsplDictLocation is the locations of the knots of each B-spline.

% Construct the level-obj.BsplScaleMin B-spline vector. The support is of 
% length 2^(obj.BsplScaleMin+1).
normalisation = 1;

LengthB = 2^(obj.BsplScaleMin-1);

numBspl = obj.sizeRT / LengthB + 3; 

u = (0 : LengthB ) / LengthB;

u = u(1:(end-1));

% B splines:
Bu3 = u.^3/6;

Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;

Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;

Bu0 = (1-u).^3/6;

% The L1 norm of the B-splines:
    N3210 = norm([Bu3 Bu2 Bu1 Bu0],normalisation);
    
    N0 = norm(Bu0,normalisation);
    
    N10 = norm([Bu1 Bu0],normalisation);
    
    N210 = norm([Bu2 Bu1 Bu0],normalisation);
    
    N3 = norm(Bu3, normalisation);
    
    N32 = norm([Bu3 Bu2] ,normalisation);
    
    N321 = norm([Bu3 Bu2 Bu1],normalisation);


BsplDict = zeros(obj.sizeRT, numBspl);

BsplDictLocation = zeros(5, numBspl);

% Boundary Splines:
BsplDict( 1 : LengthB , 1) = Bu0/N0;

    BsplDictLocation(1, 1) = 1;
        
    BsplDictLocation(2, 1) = LengthB + 1;
        
BsplDict( 1 : 2*LengthB , 2) = [Bu1 Bu0]/N10;

    BsplDictLocation(1, 2) = 1;
        
    BsplDictLocation(2, 2) = LengthB + 1;
        
    BsplDictLocation(3, 2) = 2*LengthB + 1;
        
BsplDict( 1 : 3*LengthB , 3) = [Bu2 Bu1 Bu0]/N210;

    BsplDictLocation(1, 3) = 1;
        
    BsplDictLocation(2, 3) = LengthB + 1;
        
    BsplDictLocation(3, 3) = 2*LengthB + 1;
    
    BsplDictLocation(4, 3) = 3*LengthB + 1;        
        
BsplDict( end - (3*LengthB) + 1 : end  , end-2) = [Bu3 Bu2 Bu1]/N321;
        
    BsplDictLocation(2, end-2) = (numBspl-3)*LengthB + 1;
        
    BsplDictLocation(3, end-2) = (numBspl-2)*LengthB + 1;
        
    BsplDictLocation(4, end-2) = (numBspl-1)*LengthB + 1;        
        
    BsplDictLocation(5, end-2) = (numBspl)*LengthB + 1;
    
BsplDict( end - (2*LengthB) + 1 : end  , end-1) = [Bu3 Bu2]/N32;
        
    BsplDictLocation(3, end-1) = (numBspl-2)*LengthB + 1;
        
    BsplDictLocation(4, end-1) = (numBspl-1)*LengthB + 1;        
        
    BsplDictLocation(5, end-1) = (numBspl)*LengthB + 1;
    
BsplDict( end - LengthB + 1 : end , end) = Bu3/N3;
        
    BsplDictLocation(4, end) = (numBspl-1)*LengthB + 1;        
        
    BsplDictLocation(5, end) = (numBspl)*LengthB + 1;
    
for i = 4 : (numBspl-3)
    
    BsplDict(((i-4)*LengthB+1) : (i*LengthB) , i) = [Bu3 Bu2 Bu1 Bu0]/N3210;
           
    BsplDictLocation(1, i) = (i-4)*LengthB+1;
        
    BsplDictLocation(2, i) = (i-3)*LengthB + 1;
        
    BsplDictLocation(3, i) = (i-2)*LengthB + 1;
        
    BsplDictLocation(4, i) = (i-1)*LengthB + 1;        
        
    BsplDictLocation(5, i) = i*LengthB + 1;   

end

BsplDict = sparse(BsplDict);
    
BsplInterpCoef = zeros(numBspl, obj.sizeMZ, obj.numSamples);

for i = 1 : obj.numSamples
    
    BsplInterpCoef(:,:,i) = (BsplDict'*BsplDict) \ (BsplDict' * obj.LCMS_DEFORMED(:,:,i));    
    
end
