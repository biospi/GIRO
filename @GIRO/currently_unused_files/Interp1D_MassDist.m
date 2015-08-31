%function XCMS_Interp = Interp1D_MassDist(obj, BsplInterpCoef, BsplDictRepInterp, indOrg, indNew)
function XCMS_Interp = Interp1D_MassDist(obj,XCMS, indOrg, indNew, resoluteInterp)

% This function takes in the 3D matrix XCMS and the deformation field:
%
% Interp1D_MassDist(XCMS, DeformField)
%
% XCMS: 3D sample matrix with the first dimension of
% retention time, the second dimension mass-to-charge and the third
% dimension the samples. 
% DeformField: the 2D deformation field with the first dimension samples,
% and the second dimension deformations.

%% Construct the Interpolating dictionary:
% sizeRT = size(BsplInterpCoef,1);
% 
% sizeMZ = size(BsplInterpCoef,2);
% 
% numSamples = size(BsplInterpCoef,3);

sizeRT = size(XCMS,1);

sizeMZ = size(XCMS,2);

numSamples = size(XCMS,3);

%%

XCMS_Interp = zeros(sizeRT, sizeMZ, numSamples);

XCMS_BsplRep = zeros(sizeRT*resoluteInterp, sizeMZ);

%indOrg = round(indOrg);

%linearInterpFilter = 

for i = 1 : numSamples

    % Square upsample:
    tmp = XCMS(:,:,i)/resoluteInterp;
    
    for j = 1 : resoluteInterp
     
        XCMS_BsplRep(j:resoluteInterp:end) = tmp; 

    end 
    
    % Linear upsample:
    
    
    for j = 2 : sizeRT-1 
        
        ind = (indOrg(i,:) <= indNew(j)) & (indOrg(i,:) >= indNew(j-1));
            
%        find(ind~=0)
        
        XCMS_Interp(j,:,i) = sum(XCMS_BsplRep(ind,:), 1);
        
%         if j == 2801
%             
%             s = 1;
%             
%         end
        
    end
    
    % Boundary: beginning:
    ind = ((indOrg(i,:) < indNew(1))) & (indOrg(i,:) > 0); % Throw away the warped-out

    XCMS_Interp(1,:,i) = sum(XCMS_BsplRep(ind,:),1); 
    
    if sum(ind) == 0
        
        XCMS_Interp(1,:,i) = XCMS_Interp(2,:,1);
        
    end
        
    % Boundary: end:
    ind = (indOrg(i,:) >= indNew(end-1)) & (indOrg(i,:) <= indNew(end)); % Throuw away the warped-out
    
    XCMS_Interp(end,:,i) = sum(XCMS_BsplRep(ind,:),1);
%     
    if sum(ind) == 0
        
        XCMS_Interp(end,:,i) = XCMS_Interp(end-1,:,1);
        
    end
        
end

% figure(6)
% 
% clf
% 
% plot(TICOrg')
% 

 
end




