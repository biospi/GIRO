function normalise_median(OBJ_Normalise)

N_Median = zeros(OBJ_Normalise.numSamples, 1);

for i = 1 : OBJ_Normalise.numSamples
    
    tmp = OBJ_Normalise.NMASK(:,:,i);
    N_Median(i) = median(tmp(:));
    
end

OBJ_Normalise.NMASK = N_Median;