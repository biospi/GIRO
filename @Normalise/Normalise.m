classdef Normalise
    
    properties (Access = private)
        
        Samples
                
        numSamples
        
        NMASK
          
        RefNorm % Target for normalisation
        
        TarNorm % Source for normalisation
        
    end
    
    methods (Access = private)
               
        % Generate a multi-resolution B-spline dictionary for normalisation:
%        [numBspl_1D, numBspl_2D, B] = get_BsplDict_norm(OBJ_GIRO, minScale1D, maxScale1D, minScale2D, maxScale2D, sizeMZ_Channel, normalisation);
        function BsplBasis = get_BsplDict_norm_LS(lengthB)
        
            sizeRT = size(Samples,1);
            
            sizeDyadicRT = pow2( ceil( log2( sizeRT )));
            indRT_Start = ceil( (sizeDyadicRT - sizeRT) / 2);
            
            u = (0 : lengthB ) / lengthB;
            u = u(1:(end-1));
            
            Bu3 = u.^3/6;
            Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
            Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;
            Bu0 = (1-u).^3/6;
            
            Bu = [Bu3 Bu2 Bu1 Bu0];
            Bu = Bu / norm(Bu,1); % L_inf normalise

            lengthB = length(Bu)/4;

            % With boundary:
            numBasis = sizeDyadicRT/lengthB + 3;

            % Without boundary:
            %numBasis = sizeRT/lengthB - 3;

            BsplBasis = zeros(sizeDyadicRT, numBasis);

            BsplBasis( 1 : lengthB, 1) = Bu(3*lengthB+1 : 4*lengthB);  

            BsplBasis( 1 : 2*lengthB, 2) = Bu(2*lengthB+1 : 4*lengthB);

            BsplBasis( 1 : 3*lengthB, 3) = Bu(lengthB+1 : 4*lengthB);

            BsplBasis( end - lengthB + 1 : end,end) = Bu(1 : lengthB);

            BsplBasis( end - 2*lengthB + 1 : end,end-1) = Bu(1 : 2*lengthB);

            BsplBasis( end - 3*lengthB + 1 : end, end-2) = Bu(1 : 3*lengthB);

            for k = 4 : numBasis-3

                BsplBasis((k-4) * lengthB + 1 : (k) * lengthB ,k) = Bu;
    
            end
            
            BsplBasis = sparse(BsplBasis(indRT_Start : (indRT_Start + sizeRT-1), :));
            
        end

        
        
    end
    
    methods (Access = public)
        
        function OBJ_Normalise = Normalise(OBJ_Data)
           
            OBJ_Normalise.Samples = OBJ_Data.get_Samples();
            
            OBJ_Normalise.numSamples = OBJ_Data.get_numSamples();
            
        end
            
        
        % Normalisation schemes:
        NMASK = normalise_LS(OBJ_Normalise, lengthB)
        
        N_Median = normalise_median(OBJ_Normalise)
        
        function OBJ_Normalise = update_normalise(OBJ_Normalise, Samples)

            OBJ_Normalise.Samples = Samples;
            
            OBJ_Normalise.numSamples = size(Samples,3);
                        
        end
            
        
    end
    
end

   