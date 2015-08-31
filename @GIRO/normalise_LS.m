function NMASK = normalise_LS(obj, Samples, LengthBu)
% Deriving the normalising factors:
        
LogSamples = log(Anscombe(Samples));               
 
meanLogSamples = mean(LogSamples,3);

meanLogSamples = repmat(meanLogSamples,[1,1,obj.numSamples]);

NFactor = meanLogSamples - LogSamples;

sizeRT = size(Samples,1);

%LengthBu = 2; % a quarter of B-splines support

u = (0 : LengthBu ) / LengthBu;
u = u(1:(end-1));
            
Bu3 = u.^3/6;
Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;
Bu0 = (1-u).^3/6;
            
Bu = [Bu3 Bu2 Bu1 Bu0];
Bu = Bu / norm(Bu,1); % L_inf normalise

LengthBu = length(Bu)/4;

% With boundary:
numBasis = sizeRT/LengthBu + 3;

% Without boundary:
%numBasis = sizeRT/LengthBu - 3;

BsplBasis = zeros(sizeRT, numBasis);

BsplBasis( 1 : LengthBu, 1) = Bu(3*LengthBu+1 : 4*LengthBu);  

BsplBasis( 1 : 2*LengthBu, 2) = Bu(2*LengthBu+1 : 4*LengthBu);

BsplBasis( 1 : 3*LengthBu, 3) = Bu(LengthBu+1 : 4*LengthBu);

BsplBasis( end - LengthBu + 1 : end,end) = Bu(1 : LengthBu);

BsplBasis( end - 2*LengthBu + 1 : end,end-1) = Bu(1 : 2*LengthBu);

BsplBasis( end - 3*LengthBu + 1 : end, end-2) = Bu(1 : 3*LengthBu);

for k = 4 : numBasis-3

    BsplBasis((k-4) * LengthBu + 1 : (k) * LengthBu ,k) = Bu;
    
end

BsplBasis = sparse(BsplBasis);
 
% Fit the NFactor by least square:
NMASK = zeros(size(NFactor));

for k = 1 : obj.numSamples

    LS_APP = BsplBasis \ NFactor(:,:,k);
    
    NMASK(:,:,k) = exp(2 * BsplBasis * LS_APP); % The factor of 2 compensates the Anscombe
    
end
            
end
   