% Function for generating a separate input file, to be used in a chemical
% equilibrium simulation.'
% 
% We will only use three points in x, with no-flux
% conditions at the boundaries.
% 
% A single layer is sufficient, and no transport occurs.
function chemEqInput = generateChemEqInputFile(folderName, fileName, layerNumber)

    chemFileName = 'chemEqInputFile.json';
    fullFileName = [folderName fileName];

    fileContents = fileread(fullFileName);
    chemFileOutput = jsondecode(fileContents);

    chemFileOutput.newStart = 1;
    chemFileOutput.doOutput = 0;
    chemFileOutput.doPlot = 0;
    chemFileOutput.dtPlot = 1e-11;
    
    chemFileOutput.dt = 1e2;
    chemFileOutput.tStart = 0;
    chemFileOutput.tEnd = 1e5;

    chemFileOutput.voltageSweepOn = 0;
    chemFileOutput.appliedVoltage = 0;

    chemFileOutput.interfaceLocs = [0; 1e-6];

    chemFileOutput.leftElectrodeBC = 1;
    chemFileOutput.rightElectrodeBC = 1;
    chemFileOutput.leftIonExchangeCurrent = zeros(size(chemFileOutput.leftIonExchangeCurrent));
    chemFileOutput.rightIonExchangeCurrent = zeros(size(chemFileOutput.rightIonExchangeCurrent));

    chemFileOutput.layers = chemFileOutput.layers(layerNumber);
    chemFileOutput.layers.dxMax = 2e-7;
    chemFileOutput.layers.dxMin = 2e-7;
    chemFileOutput.layers.gridSymmetry = 0;
    chemFileOutput.layers.doChemEq = 0;

    nSpeciesLayer = length(chemFileOutput.layers.species);
    for sIndex = 1:nSpeciesLayer
        chemFileOutput.layers.species(sIndex).diffusionCoeff = 0;
%         chemFileOutput.layers.species(sIndex).valence = 0;
    end

    nReactionsLayer = length(chemFileOutput.layers.reactions);
    for rIndex = 1:nReactionsLayer
        chemFileOutput.layers.reactions(rIndex).wienCoeff = 0;
    end

    outputText = jsonencode(chemFileOutput, 'PrettyPrint', true);

    fid = fopen([folderName chemFileName], 'w');
    fprintf(fid, outputText);
    fclose(fid);

    chemEqInput = chemFileName;
end