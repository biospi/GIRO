function OBJ_Normalise = normalise_LS(OBJ_Normalise, lengthB)
% Deriving the normalising factors:
        
OBJ_Normalise.TarNorm = log(GIRO.Anscombe(OBJ_Normalise.Samples));               
 
OBJ_Normalise.RefNorm = mean(OBJ_Normalise.TarNorm,3);

BsplBasis = OBJ_Normalise.get_BsplDict_norm_LS(lengthB);

for i = 1 : OBJ_Normalise.numSamples

    LS_APP = BsplBasis \ (OBJ_Normalise.RefNorm - OBJ_Normalise.TarNorm(:,:,i));
    
    OBJ_Normalise.NMASK(:,:,i) = exp(2 * BsplBasis * LS_APP); % The factor of 2 compensates the Anscombe
    
end
            
end
   