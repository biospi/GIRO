function NMASK = normalise_LS(OBJ_Normalise, lengthB)
% Deriving the normalising factors:
        
sizeRT = size(OBJ_Normalise.Samples,1);

OBJ_Normalise.TarNorm = log(Anscombe(OBJ_Normalise.Samples));               
 
OBJ_Normalise.RefNorm = mean(OBJ_Normalise.TarNorm,3);

OBJ_Normalise.RefNorm = repmat(OBJ_Normalise.RefNorm,[1,1,OBJ_Normalise.numSamples]);

NFactor = OBJ_Normalise.RefNorm - OBJ_Normalise.TarNorm;

BsplBasis = get_BsplDict_norm_LS(sizeRT, lengthB);
 
% Fit the NFactor by least square:
NMASK = zeros(size(NFactor));

for k = 1 : OBJ_Normalise.numSamples

    LS_APP = BsplBasis \ NFactor(:,:,k);
    
    NMASK(:,:,k) = exp(2 * BsplBasis * LS_APP); % The factor of 2 compensates the Anscombe
    
end
            
end
   