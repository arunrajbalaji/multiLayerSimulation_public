% genOverallMeshConstantEDLRes: Generates 1D non-uniform mesh points for a
% multi-layer geometry, by calling on appropriate individual layer mesh
% generation functions. Reads parameters that are provided as arrays, with
% an entry in each array for each layer.
%
% This version of the mesh generation function ensures that there is
% constant resolution in the EDL regions near each electrode. Basically,
% the EDL size is given, and this entire region is resolved with a
% constant mesh of size dxMin. Stretching toward dxMax is applied outside
% of the EDL region.
%
% [xCenter, xFace, dxC, dxF] = genOverallMesh(xLocs, dxMaxs, dxMins, LIDs, LRCs, EDLSize)
%
% Inputs:
%       xLocs       - x - Locations of interfaces in increasing order,
%                     including domain boundaries. Size one larger than
%                     number of layers.
%       dxMaxs      - Maximum allowable mesh spacing (array, entry per layer)
%       dxMins      - Minimum allowable mesh spacing (array, entry per layer)
%       LIDs        - Layer identification numbesr (not used currently, but
%                     provided for covenience).
%       LRCs        - Left, right, center nonuniformity (array, entry per
%                     layer)
%       EDLSize     - Size (in meters) of the EDL region for each layer.
%                     This region is resolved using a uniform mesh with
%                     size dxMin
%
% Outputs:
%       xCenter    - array of ordered mesh cell centers
%       xFace      - array of ordered mesh cell faces      
%       dxC        - cell center grid sizes
%       dxF        - cell face grid sizes
%
%
% Other m-files required: genLayerMesh.m
% MAT-files required: none
%
% See also: 
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2021
%------------- BEGIN CODE --------------
function [xCenter, xFace, dxC, dxF] = genOverallMeshConstantEDLRes(xLocs, dxMaxs, dxMins, LIDs, LRCs, EDLSize)
% Allocate vectors for output
xCenter = [];
xFace = [];

% Form input vectors for individual layer mesh generation function
xStarts = xLocs(1:end-1);
xEnds = xLocs(2:end);

for ii = 1:(length(dxMaxs))

    % Perform individual layer mesh generation
    [x_center_raw, x_face_raw] ...
        = genLayerMeshConstantEDLRes(xStarts(ii), xEnds(ii), dxMaxs(ii), dxMins(ii), LIDs(ii), LRCs(ii), EDLSize(ii));

    % Assign centers
    xCenter = [xCenter; x_center_raw];

    % Ignore last value each time, since same as first of next segment.
    xFace = [xFace; x_face_raw(1:end-1)];
end

% Add last value
xFace = [xFace; xEnds(end)];

dxC = xFace(2:end) - xFace(1:end-1);
dxF = xCenter(2:end) - xCenter(1:end-1);

%------------- END OF CODE --------------
