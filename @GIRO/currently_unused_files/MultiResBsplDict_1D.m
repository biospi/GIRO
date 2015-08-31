function [indU Dict gradDict HDict] = MultiResBsplDict_1D(minScale, maxScale, ResoluteRT,normalisation)

% This function Dict = BsplDict(minScale, maxScale, ScaleRT) calculate and
% return a dictionary of multi-resolution B splines.

ScaleBspl = minScale:maxScale; % the spacing of the coasest grid in RT dim: 2^SpacingRT

numBspl = 0;

numBsplLevelN = zeros(length(ScaleBspl),1);

for i = ScaleBspl
   
    LengthB = 2^(i-2);
    
    numBspl = numBspl + ResoluteRT / LengthB + 3; 
        
    numBsplLevelN(i) = ResoluteRT / LengthB + 3;
    
end

indU = zeros(5, numBspl);

Dict = zeros(ResoluteRT, numBspl);

gradDict = zeros(ResoluteRT, numBspl);

HDict = zeros(ResoluteRT, numBspl);
    
counter = 1;
    
for i = ScaleBspl
    
    % The underlying deformation field: initially linear
    u = (0 : 2^(i-2) ) / (2^(i-2) );
    
    u = u(1:(end-1));
    

    % B splines:
    Bu3 = u.^3/6;

    Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;

    Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;

    Bu0 = (1-u).^3/6;

    N3210 = norm([Bu3 Bu2 Bu1 Bu0],normalisation);
    
    N0 = norm(Bu0,normalisation);
    
    N10 = norm([Bu1 Bu0],normalisation);
    
    N210 = norm([Bu2 Bu1 Bu0],normalisation);
    
    N3 = norm(Bu3, normalisation);
    
    N32 = norm([Bu3 Bu2] ,normalisation);
    
    N321 = norm([Bu3 Bu2 Bu1],normalisation);
    
    % Gradient of B splines:
    dBu3 = u.^2/2;

    dBu2 = -1.5*u.^2 + u + .5;

    dBu1 = 1.5*u.^2 - 2*u;

    dBu0 = -(1-u).^2/2;
    
    % Hessian of B splines:    
    hBu3 = u;

    hBu2 = -3*u + 1;

    hBu1 = 3*u - 2;

    hBu0 = 1-u;
          
    LengthB = 2^(i-2);
    
    numBspl =  ResoluteRT / LengthB - 3; 
        
% %     % Boundary:

% Indexes of different sections of B spline
    indU(1, counter) = 1;
        
    indU(2, counter) = LengthB + 1;
        
    indU(3, counter) = 0;
        
    indU(4, counter) = 0;        
        
    indU(5, counter) = 0;     
    
    Dict( 1 : LengthB , counter) = Bu0 / N0;
 
    gradDict( 1 : LengthB , counter) = dBu0 /N0;
    
    HDict( 1 : LengthB , counter) = hBu0 / N0;
    
    counter = counter + 1;
        
    indU(1, counter) = 1;
        
    indU(2, counter) = LengthB + 1;
        
    indU(3, counter) = 2*LengthB + 1;
        
    indU(4, counter) = 0;        
        
    indU(5, counter) = 0;     
    
    Dict( 1 : 2*LengthB , counter) = [Bu1 Bu0] / N10;
        
    gradDict( 1 : 2*LengthB , counter) = [dBu1 dBu0]  / N10;
    
    HDict( 1 : 2*LengthB , counter) = [hBu1 hBu0]  / N10;
    
    counter = counter + 1;
        
    indU(1, counter) = 1;
        
    indU(2, counter) = LengthB + 1;
        
    indU(3, counter) = 2*LengthB + 1;
        
    indU(4, counter) = 3*LengthB + 1;        
        
    indU(5, counter) =    0;     
    
    Dict( 1 : 3*LengthB , counter) = [Bu2 Bu1 Bu0] / N210;
    
    gradDict( 1 : 3*LengthB , counter) = [dBu2 dBu1 dBu0] / N210;
    
    HDict( 1 : 3*LengthB , counter) = [hBu2 hBu1 hBu0] / N210;
        
    counter = counter + 1;
                
    % Middle:
    for j = 1 : numBspl  
                
        % Indexes of different sections of B spline
        indU(1, counter) = (j-1)*LengthB+1;
        
        indU(2, counter) = j*LengthB + 1;
        
        indU(3, counter) =  (j+1)*LengthB + 1;
        
        indU(4, counter) =    (j+2)*LengthB + 1;        
        
        indU(5, counter) =    (j+3)*LengthB + 1;      
        
        Dict( ((j-1)*LengthB+1) : ((j+3)*LengthB) , counter) = [Bu3 Bu2 Bu1 Bu0] / N3210;
 
        gradDict( ((j-1)*LengthB+1) : ((j+3)*LengthB) , counter) = [dBu3 dBu2 dBu1 dBu0] / N3210;
        
        HDict( ((j-1)*LengthB+1) : ((j+3)*LengthB) , counter) = [hBu3 hBu2 hBu1 hBu0] / N3210;
        
        counter = counter + 1;
     
    end
    
%     % Boundary:
    
    Dict( end - (3*LengthB) + 1 : end , counter) = [Bu3 Bu2 Bu1] / N321;
    
    % Indexes of different sections of B spline
    indU(1, counter) = 0;
        
    indU(2, counter) = j*LengthB + 1;
        
    indU(3, counter) = (j+1)*LengthB + 1;
        
    indU(4, counter) = (j+2)*LengthB + 1;        
        
    indU(5, counter) = (j+3)*LengthB + 1;      
 
    gradDict( end - (3*LengthB) + 1 : end , counter) = [dBu3 dBu2 dBu1] / N321;
        
    HDict( end - (3*LengthB) + 1 : end , counter) = [hBu3 hBu2 hBu1] / N321;
    
    counter = counter + 1;
                
    Dict( end - (2*LengthB) + 1 : end , counter) = [Bu3 Bu2] / N32;
    
    % Indexes of different sections of B spline
    indU(1, counter) = 0;
        
    indU(2, counter) = 0;
        
    indU(3, counter) =  (j+1)*LengthB + 1;
        
    indU(4, counter) =    (j+2)*LengthB + 1;        
        
    indU(5, counter) =    (j+3)*LengthB + 1;      
    
    gradDict( end - (2*LengthB) + 1 : end , counter) = [dBu3 dBu2] / N32;
        
    HDict( end - (2*LengthB) + 1 : end , counter) = [hBu3 hBu2] / N32;    
    
    counter = counter + 1;
            
    Dict( end - LengthB + 1 : end , counter) = Bu3 / N3;
    
    % Indexes of different sections of B spline
    indU(1, counter) = 0;
        
    indU(2, counter) = 0;
        
    indU(3, counter) = 0;
        
    indU(4, counter) = (j+2)*LengthB + 1;        
        
    indU(5, counter) = (j+3)*LengthB + 1;      
    
    gradDict( end - LengthB + 1 : end , counter) = dBu3 / N3;
        
    HDict( end - LengthB + 1 : end , counter) = hBu3 / N3;
           
    counter = counter + 1;
%     
end

Dict = sparse(Dict);

gradDict =  sparse(gradDict);

HDict = sparse(HDict);
