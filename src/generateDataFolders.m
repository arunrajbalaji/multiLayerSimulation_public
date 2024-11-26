function newFolderPathList = generateDataFolders(varargin)
% parseInputFile: Reads input file in JSON format, parses file, and
% generates two objects that store information. One stores general
% information about the simulation, and the other stores information about
% each layer specifically. The layer info is stored as layer objects inside
% an array.
%
% [layerInfo, uniqueSpeciesNames] = parseInputFile(filePath)
%
% Inputs:
%       filePath    - Full path to input file
%
% Outputs:
%       x_center    - array of ordered mesh cell centers
%       x_face      - array of ordered mesh cell faces      
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also:
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 16-June-2021
%------------- BEGIN CODE --------------
jobNumber = varargin{3};
varargin = {varargin{1}, varargin{2}, varargin{4:end}};
if length(varargin) >= 6

    filePath = varargin{1};
    taskIndex = varargin{2};
    
    numVarsToSweep = (length(varargin)-2)/4;
    
    fileContents = fileread(filePath);
    rawOutput = jsondecode(fileContents);
    
    parameterValues = cell(1, numVarsToSweep);
    parameterNames = cell(1, numVarsToSweep);
    
    numCases = 1;
    
    for ii = 1:numVarsToSweep
        parameterNames{ii} = varargin{2+4*(ii-1)+1};
        parameterValues{ii} = ...
            varargin{3+4*(ii-1)+1}:varargin{3+4*(ii-1)+2}:varargin{3+4*(ii-1)+3};
        numCases = numCases * length(parameterValues{ii});
    end
    
    fullParameterMatrix = combvec(parameterValues{:});
    
    [folderPath,~,~] = fileparts(filePath);
    
    newFolderPathList = cell(numCases, 1);

    for ii = 1:numCases
        newFolderName = [];
        
        for jj = 1:numVarsToSweep
            newFolderName = [newFolderName parameterNames{jj} '_' num2str(fullParameterMatrix(jj, ii)) '_'];
        end
        
        newFolderPath = [folderPath '/output/' newFolderName];
        newFolderPathList{ii} = [newFolderPath '/'];
        
        if taskIndex == 1
            if ~exist(newFolderPath, 'dir')
                mkdir(newFolderPath)
            end
            for jj = 1:numVarsToSweep
                if ~strcmp(parameterNames{jj}, 'coIonActivityCoeff')
                    eval(['rawOutput.' parameterNames{jj} ' = fullParameterMatrix(jj, ii);'])
                else
                    membraneLayers = (rawOutput.membraneLayers)';
                    for kk = membraneLayers
                        coIonIndices = (rawOutput.layers(kk).coIonIndices)';
                        for ll = coIonIndices
                            rawOutput.layers(kk).species(ll).activityCoeff = fullParameterMatrix(jj, ii);
                        end
                    end
                end
            end
            
            outputText = jsonencode(rawOutput, 'PrettyPrint', true);
            
            fid = fopen([newFolderPath '/inputFile_' num2str(jobNumber) '.json'], 'w');
            fprintf(fid, outputText);
            fclose(fid);
            
            if isfile([folderPath '/startFile.bin'])
                copyfile([folderPath '/startFile.bin'], [newFolderPath '/startFile.bin']);
            end
        end
    end
    
else
    filePath = varargin{1};
    [folderPath,~,~] = fileparts(filePath);
    newFolderPathList = {[folderPath '/']};
    
    if (taskIndex == 1)
        copyfile(filePath, [folderPath '/inputFile.json']);
    end
    
end