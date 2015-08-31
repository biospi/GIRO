function [ctn, de_dc] = interf_L1General_1stOrder_deform(OBJ_GIRO, CP)
% This function is the interface between the solver and the XCMS class. It 
% is called by the solver but uses the propertie of the  

persistent iterationCounter

persistent stepSize

if sum(CP(:)) == 0
    
    iterationCounter = 0;
    
    stepSize = nan;

else
    
    iterationCounter = iterationCounter + 1;
    
end 

numCP = numel(CP)/OBJ_GIRO.numSamples;

% Reshape CP back:      
CP = reshape( CP, numCP, OBJ_GIRO.numSamples)';

BsplDictDeform = OBJ_GIRO.BsplDictDeform;

%% Construct the integrated 3rd order B-spline basis matrix and its 
%  derivative matrix: the 3rd order B-spline basis matrix
indRT = repmat(1:OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.numSamples, 1); % Interpolating sites

if sum(sum(OBJ_GIRO.indCP(:,end-numCP+1:end))) ~= 0
    
    DeformField = (BsplDictDeform * (CP .* OBJ_GIRO.indCP(:,end-numCP+1:end))')';
    
else
 
%DeformField = (BsplDictDeform * CP')';
DeformField = 0;

end 

indKnot = indRT - DeformField; % Knots
    
d = 3; % Order of the B-splines 
    
Samples_Deformed = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples);
    
dI_dg = zeros(OBJ_GIRO.sizeDyadicRT, OBJ_GIRO.sizeMZ, OBJ_GIRO.numSamples);

SeaMassCoef = OBJ_GIRO.SeaMassCoef;

NMASK = OBJ_GIRO.NMASK; 

% Deformation using B-splines:      
for j = 1 : OBJ_GIRO.numSamples
    
    [IntBspl, DBspl] = BsplIntegral(indKnot(j,:), indRT(j,:), d, 0);

    Samples_Deformed(:,:,j) = (IntBspl * SeaMassCoef(:,:,j)) .* NMASK(:,:,j);
    
    dI_dg(:,:,j) = (DBspl * SeaMassCoef(:,:,j)) .* NMASK(:,:,j);
    
end
  
% Deforming according to CP:
      
% B-spline image representation:
           
% Criterion re-evaluation: 
[ctn, de_dI] = OBJ_GIRO.get_ctn_var_LogAnscombe_2D(Samples_Deformed, OBJ_GIRO.RTWinCTN_LevelN);  
ctn

% Chain rule to combine the three: integrate across the deformation 
% field for dg_dC:
de_dc = zeros(OBJ_GIRO.numSamples, size(CP,2)); 

if sum(sum(OBJ_GIRO.indCP(:,end-numCP+1:end))) ~= 0
    
for j = 1 : OBJ_GIRO.numSamples 
     
    % sub-gradient including the threshold:
    de_dc(j,:) = sum(de_dI(:,:,j) .* dI_dg(:,:,j) ./ (OBJ_GIRO.sizeMZ*(Samples_Deformed(:,:,j)+.375)), 2)' * BsplDictDeform;
   
    de_dc(j,OBJ_GIRO.indCP(j,(end-numCP+1:end))==0) = 0;
    
end 

end

if (iterationCounter < 2) && ((max(de_dc(:)) > 1) || (min(de_dc(:)) < .1)) && (sum(CP(:)) == 0) 

    stepSize = .375 / max(de_dc(:));
    
else if isnan(stepSize)

        stepSize = 1;
         
    end
    
end

de_dc = stepSize * reshape(de_dc', numel(de_dc), 1);
          
end  