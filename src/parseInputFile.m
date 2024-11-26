function [numericalInfo, layerInfo, constants, uniqueSpeciesNames, rxnInfo] = parseInputFile(filePath)
% parseInputFile: Reads input file in JSON format, parses file, and
% generates three objects that store information. <layerInfo> is raw
% storage of user-specified input data in a layer-by-layer format.
% <uniqueSpeciesNames> is a list of unique species names in the entire
% simulation, and <rxnInfo> stores information about reactions.
%
% [layerInfo, uniqueSpeciesNames, rxnInfo] = parseInputFile(filePath)
%
% Inputs:
%       filePath    - Character Array. Full path to input file, in
%                    pre-determined format. See sample input file for
%                    formatting requirements.
%
% Outputs:
%       numericalInfo       - Structure that contains numerical parameters
%                             that pertain to the entire simulation
%       layerInfo           - Structure Array that contains parameters and
%                             constants that vary from layer to layer.
%       constants           - Structure that contains physical constants
%                             that affect all layers. Eventually, this
%                             structure is augmented to contain physical
%                             constant arrays for those constants that
%                             vary in space and are involved in vectorized
%                             computations. See genPhysConstArray().
%       uniqueSpeciesNames  - Cell array of unique species in simulation
%       rxnInfo             - Structure array used to store reaction
%                             information for each reaction, for each layer.
%
% Other m-files required: none
% MAT-files required: MAT file may be required if electric potential
%                     boundary condition is a user-defined function of time.
%
% See also: genPhysConstArrays.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2022
%------------- BEGIN CODE --------------
% Read JSON file and store raw information
fileContents = fileread(filePath);
rawOutput = jsondecode(fileContents);

% Package numerical information pertaining to the entire simulation
numericalInfo = struct();

numericalInfo.newStart = rawOutput.newStart;
numericalInfo.restartFileName = rawOutput.restartFileName;
numericalInfo.highPrecision = rawOutput.highPrecision;
numericalInfo.numberOfDigits = rawOutput.numberOfDigits;

numericalInfo.doOutput = rawOutput.doOutput;
numericalInfo.dtOut = rawOutput.dtOut;
numericalInfo.doPlot = rawOutput.doPlot;
numericalInfo.dtPlot = rawOutput.dtPlot;

numericalInfo.dt = rawOutput.dt;
numericalInfo.tStart = rawOutput.tStart;
numericalInfo.tEnd = rawOutput.tEnd;

numericalInfo.appliedVoltage = rawOutput.appliedVoltage;
numericalInfo.voltageFileName = rawOutput.voltageFileName;
numericalInfo.voltageSweepRate = rawOutput.voltageSweepRate;
numericalInfo.voltageSweepOn = rawOutput.voltageSweepOn;

numericalInfo.interfaces = rawOutput.interfaceLocs;

numericalInfo.leftElectrodeBC = rawOutput.leftElectrodeBC;
numericalInfo.rightElectrodeBC = rawOutput.rightElectrodeBC;

numericalInfo.concentrationReductionFactor = rawOutput.concentrationReductionFactor;
numericalInfo.fixBoundaryCO2 = rawOutput.fixBoundaryCO2;
numericalInfo.co2BCConc = rawOutput.co2BCConc;

% Package constants that apply to the entire simulation.
constants = struct();
constants.T = rawOutput.T;
constants.e = rawOutput.e;
constants.kb = rawOutput.kb;
constants.stericA = rawOutput.stericA;
constants.waterPerm = rawOutput.waterPermittivity;
constants.litersPerCubicMeter = rawOutput.litersPerCubicMeter;
constants.squareCmPerSquareM = rawOutput.squareCmPerSquareM;
constants.nA = rawOutput.nA;
constants.deltaStern = rawOutput.deltaStern;
constants.deltaActivity = rawOutput.deltaActivity;

% Package properties/constants that vary from layer to layer.
layerInfo = struct([]);
nLayers = length(rawOutput.layers);
for ii = 1:nLayers
    layerInfo(ii).LID = rawOutput.layers(ii).layerID;
    layerInfo(ii).poro = rawOutput.layers(ii).porosity;
    if rawOutput.layers(ii).bruggemanOverride == 1
        layerInfo(ii).tort = rawOutput.layers(ii).porosity .^ (-(rawOutput.bruggemanConstant-1));
    else
        layerInfo(ii).tort = rawOutput.layers(ii).tortuosity;
    end
    layerInfo(ii).perm = rawOutput.layers(ii).permittivity;
    layerInfo(ii).backCharge = rawOutput.layers(ii).backCharge;
    layerInfo(ii).dxMax = rawOutput.layers(ii).dxMax;
    layerInfo(ii).dxMin = rawOutput.layers(ii).dxMin;
    layerInfo(ii).gridSymmetry = rawOutput.layers(ii).gridSymmetry;
    layerInfo(ii).doChemEq = rawOutput.layers(ii).doChemEq;
    layerInfo(ii).doMembraneEq = rawOutput.layers(ii).doMembraneEq;
end

% Determine the number and order of unique species
allSpeciesNames = [];
for ii = 1:nLayers
    nSpecies = length(rawOutput.layers(ii).species);
    for jj = 1:nSpecies
        newSpecies = string(upper(rawOutput.layers(ii).species(jj).name));
        allSpeciesNames = [allSpeciesNames; newSpecies];
    end
end

uniqueSpeciesNames = unique(allSpeciesNames);
nSpeciesTotal = length(uniqueSpeciesNames);

% Package species-specific constants for each layer. This framework ensures
% that when the order of species is changed (alphabetical, due to use of
% <unique()>), the correct value for each species is associated with that
% species' new index (in alphabetical order).
for ii = 1:nLayers
    nSpecies = length(rawOutput.layers(ii).species);

    % Logical array for determining which species are in the layer
    nameCellArray = {rawOutput.layers(ii).species(:).name};
    matchedLayers = matches(uniqueSpeciesNames,...
        nameCellArray, 'IgnoreCase', true);
    layerInfo(ii).species = matchedLayers;

    % Arrays for transferring specific constants
    % Default value zero for diffusvity (species not present)
    % Default value one for activity (species not present)
    layerInfo(ii).acti = ones(length(matchedLayers),1);
    layerInfo(ii).diff = zeros(length(matchedLayers),1);
    layerInfo(ii).vale = zeros(length(matchedLayers),1);
    layerInfo(ii).init = zeros(length(matchedLayers),1);
    layerInfo(ii).wienActivityExponent = zeros(length(matchedLayers),1);

    if ii == 1
        constants.leftIonExchangeCurrent = zeros(length(matchedLayers),1);
        constants.rightIonExchangeCurrent = zeros(length(matchedLayers),1);
    end

    for jj = 1:nSpecies
        % Determine where in the species-array we can find this specific
        % species in this layer.
        matchedLayers = matches(uniqueSpeciesNames,...
            nameCellArray{jj}, 'IgnoreCase', true);

        % Add the properties corresponding to this species in this  layer.
        layerInfo(ii).acti(matchedLayers) = [rawOutput.layers(ii).species(jj).activityCoeff]';
        layerInfo(ii).diff(matchedLayers) = [rawOutput.layers(ii).species(jj).diffusionCoeff]';
        layerInfo(ii).vale(matchedLayers) = [rawOutput.layers(ii).species(jj).valence]';
        layerInfo(ii).init(matchedLayers) = [rawOutput.layers(ii).species(jj).initVal]';
        layerInfo(ii).wienActivityExponent(matchedLayers) = [rawOutput.layers(ii).species(jj).wienActivityExponent]';
        
        if ii == 1
            constants.leftIonExchangeCurrent(matchedLayers) = [rawOutput.leftIonExchangeCurrent(jj)]';
        end
        if ii == nLayers
            constants.rightIonExchangeCurrent(matchedLayers) = [rawOutput.rightIonExchangeCurrent(jj)]';
        end
    end
end

% Fill in reaction information for each layer-for each layer: for each
% reaction in each layer, we use a 2D array to store the one or two
% reactants, one or two products, rate coefficient, and coefficient for the
% second Wien effect (zero for no second Wien effect) for each reaction.
% Repeated reactions are added to the structure all the same, no checking
% for that is completed. See genPhysConstArrays for handling of multiple
% identical reactions in same layer.
rxnInfo = struct([]);
for ii = 1:nLayers
    nReacts = length(rawOutput.layers(ii).reactions);
    rxnInfo(ii).speciesRxns = zeros(nReacts, 6);
    for jj = 1:nReacts
        reactantIndex = find(matches(uniqueSpeciesNames,...
            rawOutput.layers(ii).reactions(jj).reactants{1}, 'IgnoreCase', true) == 1);
        if ~isempty(reactantIndex)
            rxnInfo(ii).speciesRxns(jj,1) = reactantIndex;
        else
            rxnInfo(ii).speciesRxns(jj,1) = 0;
        end

        if length(rawOutput.layers(ii).reactions(jj).reactants) < 2
            rxnInfo(ii).speciesRxns(jj,2) = 0;
        else
            reactantIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).reactants{2}, 'IgnoreCase', true) == 1);
            if ~isempty(reactantIndex)
                rxnInfo(ii).speciesRxns(jj,2) = reactantIndex;
            else
                rxnInfo(ii).speciesRxns(jj,2) = 0;
            end
        end

        productIndex = find(matches(uniqueSpeciesNames,...
            rawOutput.layers(ii).reactions(jj).products{1}, 'IgnoreCase', true) == 1);
        if ~isempty(productIndex)
            rxnInfo(ii).speciesRxns(jj,3) = productIndex;
        else
            rxnInfo(ii).speciesRxns(jj,3) = 0;
        end

        if length(rawOutput.layers(ii).reactions(jj).products) < 2
            rxnInfo(ii).speciesRxns(jj,4) = 0;
        else
            productIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).products{2}, 'IgnoreCase', true) == 1);
            if ~isempty(productIndex)
                rxnInfo(ii).speciesRxns(jj,4) = productIndex;
            else
                rxnInfo(ii).speciesRxns(jj,4) = 0;
            end
        end

        rxnInfo(ii).speciesRxns(jj,5) = rawOutput.layers(ii).reactions(jj).forwardRateCoeff;
        rxnInfo(ii).speciesRxns(jj,6) = rawOutput.layers(ii).reactions(jj).wienCoeff;
    end
end

