function produceOutput(xCenter, uniqueSpecies, z, t, folderName, isHighPrecision)
% doOutput: Output simulation data at a certain time step, as a binary
% file.
%
% doOutput(xCenter, xFace, uniqueSpecies, z, t, folderName)
%
% Inputs:
%       uniqueSpecies   - ordered list of unique species names, from
%                         parsing function.
%       layerInfo       - layer info struct, as output by parsing function.
%       xCenter         - List of cell centers, where values are stored.
%                         Should not include an interface between layers.
%
% Outputs:
%       constants.poro            - 1D array of porosity
%       constants.tort            - 1D array of tortuosity
%       constants.perm            - 1D array of electric permittivity
%       constants.diff            - 2D array of species diffusivities
%       constants.acti            - 2D array of species tortuosities
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also: parseInputFile.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 25-June-2021
%------------- BEGIN CODE --------------

% Create file for writing
if ~isHighPrecision
    fileName = [folderName 'time_' num2str(t) '.bin'];
    fileID = fopen(fileName, 'w');
end

% Create data structure for writing
nSpecies = length(uniqueSpecies);
if isHighPrecision
    dataOutputMatrix = mp(zeros(nSpecies+2, length(xCenter)));
else
    dataOutputMatrix = zeros(nSpecies+2, length(xCenter));
end

dataOutputMatrix(1,:) = xCenter';

for ii = 1:nSpecies
    dataOutputMatrix(1+ii,:) = z(ii:(1+nSpecies):end)';
end

dataOutputMatrix(nSpecies+2,:) = z((nSpecies+1):(1+nSpecies):end)';

% Write data
if isHighPrecision
    mp.write(dataOutputMatrix, [folderName 'time_' num2str(t) '.txt'])
else
    fwrite(fileID, dataOutputMatrix, 'double');
    
    % Close file
    fclose(fileID);
end

%------------- END OF CODE --------------
