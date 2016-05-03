function [ctn, de_dc] = interf_register_solver(OBJ_GIRO, CP)
% This function is the interface between the solver and the register in GIRO class. 
Down_Samples = cell2mat(OBJ_GIRO.Down_Samples(OBJ_GIRO.currentLevel - OBJ_GIRO.ResRT_Begin + 1));

numCP = numel(CP)/OBJ_GIRO.numSamples;

% Reshape CP back:      
CP = reshape( CP, numCP, OBJ_GIRO.numSamples)';

%% Construct the integrated 3rd order B-spline basis matrix and its 
%  derivative matrix: the 3rd order B-spline basis matrix
indRT = repmat(1:OBJ_GIRO.sizeDyadicRT_currentLevel, OBJ_GIRO.numSamples, 1); % Interpolating sites

if sum(sum(OBJ_GIRO.indCP(:,end-numCP+1:end))) ~= 0
    
    DeformField = (OBJ_GIRO.BsplDictDeform * (CP .* OBJ_GIRO.indCP(:,end-numCP+1:end))')';
    
else
 
%DeformField = (BsplDictDeform * CP')';
DeformField = 0;

end 

indKnot = indRT + DeformField; % Knots
    
d = 3; % Order of the B-splines 
        
dI_dg = zeros(size(OBJ_GIRO.DeformedSamples));

%[GridX, GridY] = meshgrid(1:OBJ_GIRO.sizeMZ, 1:OBJ_GIRO.sizeDyadicRT);

% Deformation using B-splines:      
for j = 1 : OBJ_GIRO.numSamples

%    SeaMassCoef = filter([.0417 .4583 .4583 .0417], 1, OBJ_GIRO.SeaMassCoef(:,:,j), [], 1);
    
    [IntBspl, DBspl] = OBJ_GIRO.BsplIntegral(indKnot(j,:), indRT(j,:), d, 0);
%    [indX_SeaMassCoef, indY_SemaMassCoef] = meshgrid(1:OBJ_GIRO.sizeMZ, indKnot(j,:));

    OBJ_GIRO.DeformedSamples(:,:,j) = (IntBspl * Down_Samples(:,:,j));
%    OBJ_GIRO.DeformedSamples(:,:,j) = interp2(indX_SeaMassCoef, indY_SemaMassCoef, SeaMassCoef, GridX, GridY);

    dI_dg(:,:,j) = DBspl * Down_Samples(:,:,j);
%    dI_dg(1:end-1,:,j) = diff(OBJ_GIRO.DeformedSamples(:,:,j), 1);

end

NMASK = OBJ_GIRO.OBJ_Normalise.get_Normaliser();

OBJ_GIRO.DeformedSamples = OBJ_GIRO.DeformedSamples .* NMASK;
    
dI_dg = dI_dg .* NMASK;

% Deforming according to CP:
      
% B-spline image representation:
           
% Criterion re-evaluation: 
[ctn, de_dI] = OBJ_GIRO.OBJ_CTN_LS_2D.get_Criterion(OBJ_GIRO.DeformedSamples);  
ctn

% Chain rule to combine the three: integrate across the deformation 
% field for dg_dC:
de_dc = zeros(OBJ_GIRO.numSamples, size(CP,2)); 

if sum(sum(OBJ_GIRO.indCP(:,end-numCP+1:end))) ~= 0
    
for j = 1 : OBJ_GIRO.numSamples 
     
    % sub-gradient including the threshold:
    de_dc(j,:) = sum(de_dI(:,:,j) .* dI_dg(:,:,j) ./ (OBJ_GIRO.sizeMZ*(OBJ_GIRO.DeformedSamples(:,:,j)+.375)), 2)' * OBJ_GIRO.BsplDictDeform(OBJ_GIRO.indRT_Start_currentLevel : (OBJ_GIRO.sizeRT_currentLevel+OBJ_GIRO.indRT_Start_currentLevel-1),:);
%    de_dc(j,:) = sum(de_dI(:,:,j) .* dI_dg(:,:,j) ./ (OBJ_GIRO.DeformedSamples(:,:,j)+.375), 2)' * OBJ_GIRO.BsplDictDeform(OBJ_GIRO.indRT_Start_currentLevel : (OBJ_GIRO.sizeRT_currentLevel+OBJ_GIRO.indRT_Start_currentLevel-1),:);
   
    de_dc(j,OBJ_GIRO.indCP(j,(end-numCP+1:end))==0) = 0;
        
end 

end

de_dc = -OBJ_GIRO.searchingStepSize * reshape(de_dc', numel(de_dc), 1);
          
end  