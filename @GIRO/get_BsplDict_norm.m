function [numBspl_1D, numBspl_2D, B] = get_BsplDict_norm(OBJ_GIRO, minScale1D, maxScale1D, minScale2D, maxScale2D, sizeMZ_Channel, normalisation)

% This function Dict = BsplDict(minScale, maxScale, ScaleRT) calculate and
% return a dictionary of multi-resolution B splines.

% 1D basis:

ScaleBspl_1D = minScale1D : maxScale1D; % the spacing of the coasest grid in RT dim: 2^SpacingRT

ScaleBspl_2D = minScale2D : maxScale2D;

numBspl_1D = zeros(length(ScaleBspl_1D),1);

numBspl_2D = zeros(length(ScaleBspl_2D),1);

numElement_1D = zeros(length(ScaleBspl_1D),1);

numElement_2D = zeros(length(ScaleBspl_2D),1);

for i = ScaleBspl_1D
    
    LengthB_1D = 2^(i-2);
       
    numBspl_1D(i - minScale1D + 1) = ((OBJ_GIRO.sizeDyadicRT / LengthB_1D) + 3) * sizeMZ_Channel; 
        
    numElement_1D(i - minScale1D + 1) = ((OBJ_GIRO.sizeDyadicRT / LengthB_1D) + 3) * LengthB_1D * 4 * sizeMZ_Channel;
    
end 

for i = ScaleBspl_2D
    
    LengthB_2D = 2^(i-2);
    
    numBspl_2D(i - minScale2D + 1) = (OBJ_GIRO.sizeDyadicRT / LengthB_2D) + 3;
    
    numElement_2D(i - minScale2D + 1) = ((OBJ_GIRO.sizeDyadicRT / LengthB_2D) + 3) * LengthB_2D * 4 * sizeMZ_Channel;
    
end

B = spalloc(OBJ_GIRO.sizeDyadicRT*sizeMZ_Channel, sum(numBspl_1D) + sum(numBspl_2D), sum(numElement_1D) + sum(numElement_2D));
    
counter = 1;
    

for i = ScaleBspl_1D
    
for k = 1 : sizeMZ_Channel
    
    % The underlying deformation field: initially linear
    u = (0 : 2^(i-2) ) / (2^(i-2) );
    
    u = u(1:(end-1));

    % B splines:
    Bu3 = u.^3/6;

    Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;

    Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;

    Bu0 = (1-u).^3/6;

    % Norm used to normalise the B-splines:
    N = norm([Bu3 Bu2 Bu1 Bu0],normalisation);
        
    LengthB = 2^(i-2);
    
    % Boundary:
    B( (1 : LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = Bu0 / N;
     
    counter = counter + 1;
    
    B( (1 : 2*LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu1 Bu0] / N;
        
    counter = counter + 1;
    
    B( (1 : 3*LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu2 Bu1 Bu0] / N;
    
    counter = counter + 1;
                
    % Middle:
    for j = 1 : (numBspl_1D(i - minScale1D + 1)/ sizeMZ_Channel - 6) 
        
        B( (((j-1)*LengthB+1) : ((j+3)*LengthB))  + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2 Bu1 Bu0] / N;
 
        counter = counter + 1;
     
    end
    
%     % Boundary:
    
    B( (((numBspl_1D(i - minScale1D + 1)/sizeMZ_Channel - 6)  * LengthB + 1) : ((numBspl_1D(i - minScale1D + 1)/ sizeMZ_Channel - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2 Bu1] / N;
    
    counter = counter + 1;
                
    B( (((numBspl_1D(i - minScale1D + 1)/sizeMZ_Channel - 5)  * LengthB + 1) : ((numBspl_1D(i - minScale1D + 1)/ sizeMZ_Channel - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2] / N;
    
    counter = counter + 1;
            
    B( (((numBspl_1D(i - minScale1D + 1)/sizeMZ_Channel - 4)  * LengthB + 1) : ((numBspl_1D(i - minScale1D + 1)/ sizeMZ_Channel - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = Bu3 / N;
           
    counter = counter + 1;
%     
end

end

for i = ScaleBspl_2D 
    
    % The underlying deformation field: initially linear
    u = (0 : 2^(i-2) ) / (2^(i-2) );
    
    u = u(1:(end-1));

    % B splines:
    Bu3 = u.^3/6;

    Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;

    Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;

    Bu0 = (1-u).^3/6;

    % Norm used to normalise the B-splines:
%    N = norm([Bu3 Bu2 Bu1 Bu0],normalisation) * (OBJ_GIRO.indChannel(indChannel + 1) - OBJ_GIRO.indChannel(indChannel));
    N = norm([Bu3 Bu2 Bu1 Bu0],normalisation);% * sizeMZ_Channel * (OBJ_GIRO.numSamples-1);


    LengthB = 2^(i-2);
    
    for k = 1 : sizeMZ_Channel
    
    % Boundary:
    B( (1 : LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = Bu0 / N;
     
    end
    
    counter = counter + 1;
    
    for k = 1 : sizeMZ_Channel
    
    B( (1 : 2*LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu1 Bu0] / N;
        
    end
    
    counter = counter + 1;
    
    for k = 1 : sizeMZ_Channel
    
    B( (1 : 3*LengthB) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu2 Bu1 Bu0] / N;
    
    end
    
    counter = counter + 1;
                
    % Middle:
    for j = 1 : (numBspl_2D(i - minScale2D + 1) - 6)  
       
        for k = 1 : sizeMZ_Channel
        
        B( (((j-1)*LengthB+1) : ((j+3)*LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2 Bu1 Bu0] / N;
 
        end
        
        counter = counter + 1;
     
    end
    
%     % Boundary:
    
    for k = 1 : sizeMZ_Channel

    B( (((numBspl_2D(i - minScale2D + 1) - 6)  * LengthB + 1) : ((numBspl_2D(i - minScale2D + 1) - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2 Bu1] / N;
    
    end
    
    counter = counter + 1;
    
    for k = 1 : sizeMZ_Channel
    
    B( (((numBspl_2D(i - minScale2D + 1) - 5)  * LengthB + 1) : ((numBspl_2D(i - minScale2D + 1) - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = [Bu3 Bu2] / N;
    
    end
    
    counter = counter + 1;
            
    for k = 1 : sizeMZ_Channel
    
    B( (((numBspl_2D(i - minScale2D + 1) - 4)  * LengthB + 1) : ((numBspl_2D(i - minScale2D + 1) - 3) * LengthB)) + (k-1)*OBJ_GIRO.sizeDyadicRT , counter) = Bu3 / N;
           
    end
    
    counter = counter + 1;
%     
end