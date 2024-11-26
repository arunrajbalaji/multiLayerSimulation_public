% genLayerMesh: Generates 1D non-uniform mesh points for a single layer.
% This version of the mesh generation function ensures that there is
% constant resolution in the EDL regions. Outside of this region, cell
% centers are distributed non-unformly according to a hyperbolic-tan.
% stretching function:
%
%       x_f = (tanh((eta - 0.5)/meshParam) - tanh(-0.5/meshParam))
%               /(2*tanh(0.5/meshparam))*(xEnd - xStart) + xStart
%
% Here, eta is taken to be a vector of uniformly distributed points between
% 0 and 1 (including both values), of length N+1 total. x_f is the cell
% face locations, and x_c will be computed as midpoints of the faces.
%
% The LRC flag allows users to select whether the layer will have a
% central/symmetric hyperbolic tangent distribution of points, or a
% left/right sided half-hyperbolic-tangent distribution of points (which
% might be useful for domains that are adjacent to a diffusion layer).
%
% The number of grid points and the meshParam parameter are computed for
% the hyperbolic tangent mapping function automatically, based on the user
% specified dxMin and dxMax.
%
% This version of the mesh generation function ensures that there is
% constant resolution in the EDL regions near each electrode. Basically,
% the EDL size is given, and this entire region is resolved with a
% constant mesh of size dxMin. Stretching toward dxMax is applied outside
% of the EDL region.
%
% Syntax:  [xCenter, xFace] = genLayerMesh(xStart, xEnd, dxMax, dxMin, LID, LRC)
%
% Inputs:
%       xStart      - location of first cell face
%       xEnd        - location of last cell face
%       dxMax       - Maximum cell size, same units as xStart and xEnd
%       dxMin       - Minimum cell size, same units as xStart and xEnd
%       LID         - Layer ID Number
%       LRC         - Left sided (1), right sided (2), or central (0)
%                     hyperbolic mapping. Side refers to that with the
%                     concentrated nonuniform points.
%       EDL         - Size (in meters) of the EDL region for each layer.
%                     This region is resolved using a uniform mesh with
%                     size dxMin
%
% Outputs:
%       xCenter     - array of ordered mesh cell centers, column vector
%       xFace       - array of ordered mesh cell faces, column vectorNs(
%
%
% Other m-files required: none
% MAT-files required: none
%
% See also: genOverallMesh.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 16-June-2021
%------------- BEGIN CODE --------------
function [xCenter, xFace] = genLayerMeshConstantEDLRes(xStart, xEnd, dxMax, dxMin, LID, LRC, EDL)
try 
    % Determine cell face locations by mapping uniformly spaced points
    if LRC == 0         % CENTRAL
        % Exclude the EDL zones from the stretching region.
        xStartStretching = xStart + EDL;
        xEndStretching = xEnd - EDL;

        % Compute the number of grid points and the meshParam for this layer
        epsFunc = @(z) dxMin/(xEndStretching - xStartStretching) - (tanh((z*atanh(2*tanh(0.5/z)*dxMax/(xEndStretching - xStartStretching)) - 0.5)/(z))...
            - tanh(-0.5/z))/(2*tanh(0.5/z));
        meshParam = fzero(epsFunc, 0.1);
        N = ceil(1/(meshParam*atanh(2*tanh(0.5/meshParam)*dxMax/(xEndStretching - xStartStretching))));
        
        % Uniformly spaced points, for mapping purposes
        eta = (linspace(0,1,N+1))';
        
        xFaceStretched = (tanh((eta - 0.5)/meshParam) - tanh(-0.5/meshParam))...
            /(2*tanh(0.5/meshParam))*(xEndStretching - xStartStretching) + xStartStretching;

        xFace = [(xStart:dxMin:xStartStretching-dxMin)'; xFaceStretched; (xEndStretching+dxMin:dxMin:xEnd)'];
    
    elseif LRC == 1     % LEFT SIDED
        % Exclude the EDL zones from the stretching region.
        xStartStretching = xStart + EDL;

        % Compute the number of grid points and the meshParam for this layer
        epsFunc = @(z) dxMax/(xEnd - xStartStretching) - tanh(1/z*(1 - z*atanh((1-dxMin/(xEnd - xStartStretching))*tanh(1/z))))/tanh(1/z);
        meshParam = fzero(epsFunc, 0.1);
        N = ceil((1 - meshParam*atanh((1-dxMin/(xEnd - xStartStretching))*tanh(1/meshParam)))^(-1));
        
        % Uniformly spaced points, for mapping purposes
        eta = (linspace(0,1,N+1))';
        
        xFaceStretched = flipud(-((tanh((eta)/meshParam))...
            /(tanh(1/meshParam))-1)*(xEnd - xStartStretching) + xStartStretching);

        xFace = [(xStart:dxMin:xStartStretching-dxMin)'; xFaceStretched];
    
    elseif LRC == 2     % RIGHT SIDED
        % Exclude the EDL zones from the stretching region.
        xEndStretching = xEnd - EDL;

        % Compute the number of grid points and the meshParam for this layer
        epsFunc = @(z) dxMax/(xEndStretching - xStart) - tanh(1/z*(1 - z*atanh((1-dxMin/(xEndStretching - xStart))*tanh(1/z))))/tanh(1/z);
        meshParam = fzero(epsFunc, 0.1);
        N = ceil((1 - meshParam*atanh((1-dxMin/(xEndStretching - xStart))*tanh(1/meshParam)))^(-1));
        
        % Uniformly spaced points, for mapping purposes
        eta = (linspace(0,1,N+1))';
        
        xFaceStretched = (tanh((eta)/meshParam))...
            /(tanh(1/meshParam))*(xEndStretching - xStart) + xStart;

        xFace = [xFaceStretched; (xEndStretching+dxMin:dxMin:xEnd)'];
    end
    
    % Determine cell center locations by finding midpoints of faces
    xCenter = 0.5*(xFace(2:end) + xFace(1:end-1));
    
catch
    warning(['Error in genLayerMesh, check Layer # ' num2str(LID)])
    warning('Assigning a uniform mesh for this layer')
    
    % Uniformly spaced points, for mapping purposes
    N = 100;
    eta = (linspace(0,1,N+1))';
    
    % Determine cell face locations by mapping uniformly spaced points
    xFace = eta*(xEnd - xStart) + xStart;
    
    % Determine cell center locations by finding midpoints of faces
    xCenter = 0.5*(xFace(2:end) + xFace(1:end-1));

end

%------------- END OF CODE --------------