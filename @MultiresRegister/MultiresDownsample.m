function MultiresLCMS = MultiresDownsample(OBJ_GIRO, SeaMassCoefDyadicRT)

reduceLevels = OBJ_GIRO.SeaMassResRT_End - OBJ_GIRO.SeaMassResRT_Begin + 1;

BsplSmoothing = [1 4 1]/6;

if reduceLevels <= 0
    
    return
    
else
    
    originalImage = SeaMassCoefDyadicRT; 
    
    MultiresLCMS = cell(reduceLevels, 1);
        
    for i = 1 : reduceLevels
                
        MultiresLCMS(end - i + 1) = {originalImage}; 
        
        downsampledImage = zeros( size(originalImage, 1)/2, size(originalImage,2), OBJ_GIRO.numSamples);

        for j = 1 : OBJ_GIRO.numSamples
        
            % Smoothing using Bspl filter:
            originalImage(:,:,j) = filter(BsplSmoothing, 1, originalImage(:,:,j), [], 1);
        
            % Downsampling:
            downsampledImage(:,:,j) = originalImage(1:2:end, :, j);
                        
        end
        
        originalImage = downsampledImage;

            
    end

end