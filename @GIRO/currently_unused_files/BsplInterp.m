% function BsplDictDeform = BsplInterp(obj, BsplDeformField)

% This function receives the update of deformation field and interpolate
% the LC-MS images accordingly to generate a set of deformed B-spline bases
% for later reconstruction of LC-MS samples.
% 
% BsplDictDeform = BsplInterp(obj, BsplDeformField)
%

% Set up a look-up table for the B-spline basis:
resLUT = 32;

BsplDictLocationUpdate = obj.BsplDictLocation;

% B spline bases of resolution 4*resLUT:
u = (0 : (resLUT - 1)) / resLUT;

Bu3 = u.^3/6;

Bu2 = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;

Bu1 = ( 3*u.^3 - 6*u.^2 + 4)/6;

Bu0 = (1-u).^3/6;

B = [Bu3 Bu2 Bu1 Bu0];

% Boundary not deformed:
BoundaryWidth = 8;

for i = 1 : obj.numSamples

    % Keep a tem version of deformation field for query:
    DeformFieldI = BsplDeformField(i,:);
    
    % Deform the dictionary by querying the deform field:
    BsplDictLocationUpdate(:, BoundaryWidth:end-BoundaryWidth) = obj.BsplDictLocation(:,BoundaryWidth:end-BoundaryWidth) - DeformFieldI(obj.BsplDictLocation(:,BoundaryWidth:end-BoundaryWidth));

    % Get the integer indexes for the deformed dictionary:
    
    
    
    % L1 normaliser:
        


end

% end