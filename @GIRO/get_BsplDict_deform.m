function [BsplDictDeform, normBuL2, numCoefs] = get_BsplDict_deform(OBJ_GIRO, numResLevels)

numCoefs = zeros(1,numResLevels);
 
numNonZeroEleDeform = 0;

% Set # of total non-zero coefs and # of total non-zero elements in both
% dictionaries:
for i = 1 : numResLevels
    
    lengthB = OBJ_GIRO.lengthB_Deform * 2^(i-1);
    
    %without boundary conditions:
%    numCoefs(i) = OBJ_GIRO.sizeDyadicRT/lengthB - 3;
    
    %with boundary conditions:
    numCoefs(i) = OBJ_GIRO.sizeDyadicRT_currentLevel/lengthB + 3;
    
    numNonZeroEleDeform = numNonZeroEleDeform + numCoefs(i) * 4 * lengthB;
    
end

indLevels = [1 cumsum(numCoefs)+1];

% Allocate memory for sparse deforma basis:
BsplDictDeform = zeros(OBJ_GIRO.sizeDyadicRT_currentLevel, sum(numCoefs)); 
                                
normBuL2 = zeros(1,numResLevels);

% Deform basis:
for i = 1 : numResLevels 

    lengthB = OBJ_GIRO.lengthB_Deform * 2^(i-1);    
    
    % Construct Level i sampled B splines:
    u = .5/lengthB : 1/lengthB : (1-.5/lengthB);

%    length(u)
    
    Bu3 = u.^3/6;
    Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
    Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;
    Bu0 = (1-u).^3/6;

    Bu = [Bu3 Bu2 Bu1 Bu0]; % Look-up table for B-splines
    
    normBuL2(i) = norm(Bu,2);
    
    counterBu = 0;
    
    Bu = Bu / normBuL2(i);

    BsplDictDeform(1 : lengthB,indLevels(i)) = Bu0 / normBuL2(i);
    
    BsplDictDeform(1 : 2*lengthB,indLevels(i)+1) = [Bu1 Bu0] / normBuL2(i);
    
    BsplDictDeform(1 : 3*lengthB,indLevels(i)+2) = [Bu2 Bu1 Bu0] / normBuL2(i);
    
    for k = indLevels(i)+3 : indLevels(i+1)-4
%     for k = indLevels(i) : indLevels(i+1)-1
        
        BsplDictDeform(counterBu * lengthB + 1 : counterBu * lengthB + 4*lengthB,k) = Bu;
        
        counterBu = counterBu + 1;
    
    end
    
    BsplDictDeform(end - 3*lengthB + 1 : end,indLevels(i+1)-3) = [Bu3 Bu2 Bu1] / normBuL2(i);
        
    BsplDictDeform(end - 2*lengthB + 1 : end,indLevels(i+1)-2) = [Bu3 Bu2] / normBuL2(i);
    
    BsplDictDeform(end - lengthB + 1 : end,indLevels(i+1)-1) = Bu3 / normBuL2(i);

end

BsplDictDeform = sparse(BsplDictDeform);
          
end

