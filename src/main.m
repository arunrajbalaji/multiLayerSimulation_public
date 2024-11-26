 % main: This is the main simulation function. This function can be called
% by a script that runs several cases in parallel, if desired.
%
% main(folderName, fileName)
%
% Inputs:
%       folderName      - Character Array. Directory in which input file
%                         can be found and to which data should be saved.
%
%       fileName        - Character Array. Name of input file. See sample
%                         input file for formatting requirements.
%
% Outputs:
%       Function itself does not return outputs, data is saved as binary
%       files in the <folderName> directory if flag is turned on. See
%       <produceOutput> function for details about data output. Will also
%       produce MATLAB figures if flag is turned on in input file. See
%       <producePlot> function for details about plot production.
%
% Example: 
%       Run two cases in parallel. Each input file contained in
%       <inputFileNameArray> must be located inside its respective
%       directory in <folderArray>.
%
%       folderArray = {'/path/to/folder/one/', '/path/to/folder/two/'};
%       inputFileNameArray = {'inputFileOne.json', 'inputFileTwo.json'};
%       parfor ii = 1:length(folderArray)
%           main(folderArray{ii}, inputFileNameArray{ii})
%       end
%
% Other m-files required:
%       - parseInputFile.m
%       - genOverallMesh.m
%       - genLayerMesh.m
%       - genPhysConstArrays.m
%       - initDomain.m
%       - produceOutput.m
%       - producePlot.m
%
%       *** Advanpix Multi-precision Toolbox is required if more than
%       standard double precision is required. ***
% 
% MAT-files required:
%       MAT file may be required if electric potential boundary condition
%       is a user-defined function of time.
%
% See also: sampleInput.json
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2022
%
%------------- BEGIN CODE --------------
function main(folderName, fileName)

% Input folder and file name
fullName = [folderName, fileName]

% Parse input files and package data
[numericalInfo, layerInfo, constants, uniqueSpecies, rxnInfo] = parseInputFile(fullName);
nSpecies = length(uniqueSpecies);

% May need to perform chemical equilbiration for some layers. That is done
% here, by using a reference input file and calling the simulation from
% within itself. No transport, only chemical equilibrium. Three grid
% points, with no-flux boundary conditions.
for ii = 1:length(layerInfo)
    if layerInfo(ii).doChemEq
        chemEqInput = generateChemEqInputFile(folderName, fileName, ii);
        layerInfo(ii).init = runChemEqSimulation(folderName, chemEqInput);
    end
    % Membrane equilibration: will set concentration of solution in
    % membrane to precisely the value required for net charge neutrality,
    % considering background charge and porosity of the membrane.
    if layerInfo(ii).doMembraneEq
        memCharge = sum(layerInfo(ii).init .* layerInfo(ii).vale);
        if sign(memCharge) == sign(layerInfo(ii).backCharge)
            disp(['Error: Charge of solution in membrane ' num2str(ii) ...
                ' does not allow for charge neutrality.'])
        else
            layerInfo(ii).init = -layerInfo(ii).init * (layerInfo(ii).backCharge/layerInfo(ii).poro)/memCharge;
        end
    end
end

% Restore boundary co2 values, based on Henry's Law and continual dissolution
% into the working fluid.
co2Index = find(strcmp(uniqueSpecies, 'CO2'));
if numericalInfo.fixBoundaryCO2
    layerInfo(1).init(co2Index) = numericalInfo.co2BCConc;
    layerInfo(end).init(co2Index) = numericalInfo.co2BCConc;

    disp([uniqueSpecies{co2Index} ' concentration set using override.']);
end

% Set precision using Advanpix package if high precision required
if numericalInfo.highPrecision
    addpath('../advanpix/');
    mp.Digits(numericalInfo.numberOfDigits);
    format longG;
end

% Generate mesh for entire domain
% Last argument is a vector of double layer thickness for each layer.
% Implemented in an ad-hoc way just for testing purposes. Can be computed
% automatically based on input file. IMPROVE LATER.
% [xCenter, xFace, dxC, dxF] = genOverallMeshConstantEDLRes(numericalInfo.interfaces,...
%     [layerInfo(:).dxMax], [layerInfo(:).dxMin],...
%     [layerInfo(:).LID], [layerInfo(:).gridSymmetry], [0, 0, 0, 0, 0]);
[xCenter, xFace, dxC, dxF] = genOverallMesh(numericalInfo.interfaces,...
    [layerInfo(:).dxMax], [layerInfo(:).dxMin],...
    [layerInfo(:).LID], [layerInfo(:).gridSymmetry]);

% Physical constant array generation
constants = genPhysConstArrays(uniqueSpecies, layerInfo, constants, rxnInfo, xCenter, numericalInfo.interfaces);

% Time paramters, use user-specified values
tEnd = numericalInfo.tEnd;
dt = numericalInfo.dt;
dtOut = numericalInfo.dtOut;
dtPlot = numericalInfo.dtPlot;
tStart = numericalInfo.tStart;

totalSteps = ceil((tEnd - 0)/dt);
t = tStart;

% Save data structures that will be useful for post-processing functions.
% Write as matlab files, since this only done once - speed not an issue.
% Might be redundant, could improve later if necessary
save(strcat(folderName, '/', 'meshInfo.mat'), 'xCenter', 'xFace', 'dxC', 'dxF')
save(strcat(folderName, '/', 'layerInfo.mat'), 'layerInfo')
save(strcat(folderName, '/', 'constants.mat'), 'constants')
save(strcat(folderName, '/', 'rxnInfo.mat'), 'rxnInfo')
save(strcat(folderName, '/', 'numericalInfo.mat'), 'numericalInfo')
save(strcat(folderName, '/', 'uniqueSpecies.mat'), 'uniqueSpecies')

% Variable initialization function. Can be swapped with whatever is
% desired. Z is a vector of size [(nSpecies+1)*length(dxC), 1]. Variables
% are interwoven for each gridpoint. For species A, B, C starting at
% gridpoint 1 and ending at gridpoint end:
%
% [ A_1
%   B_1
%   C_1
%   Potential_1
%   A_2
%   B_2
%   C_2
%   Potential_2
%    .
%    .
%    .
%   A_end
%   B_end
%   C_end
%   Potential_end ]
%
% Note: gridpoints 1 and end represent ghost cells. Dirichlet boundary
% conditions are applied at ghost cell centers, flux-based boundary
% conditions are applied at the face between ghost cells and the domain.
if numericalInfo.newStart
    % Default initalization, using values from user input
    z = initDomain(nSpecies, xCenter, numericalInfo, constants);
    z(1:nSpecies) = z(1:nSpecies) * numericalInfo.concentrationReductionFactor;
    z(end-nSpecies:end-1) = z(end-nSpecies:end-1) * numericalInfo.concentrationReductionFactor;
else
    fullRestartFileName = strcat(folderName, '/', numericalInfo.restartFileName);
    if numericalInfo.highPrecision
        rawData = mp.read(fullRestartFileName);
    else
        fileID = fopen(fullRestartFileName);
        rawData = fread(fileID, 'double');
        fclose (fileID);
        nCell = length(rawData)/(nSpecies+2);
        rawData = reshape(rawData, [nSpecies+2, nCell]);
    end

    if numericalInfo.highPrecision
        z = mp(zeros((nSpecies+1)*length(xCenter), 1));
    else
        z = zeros((nSpecies+1)*length(xCenter), 1);
    end

    xCenter_old = rawData(1,:);
    for ii = 1:nSpecies
        z(ii:(nSpecies+1):end) = interp1(xCenter_old, rawData(ii+1,:), xCenter);
        z(1:nSpecies) = constants.initVal(1,:);
        z((end-nSpecies):1:end-1) = constants.initVal(end,:);
    end
%    z(nSpecies+1:nSpecies+1:end) = interp1(xCenter_old, rawData(nSpecies+2,:), xCenter);
%     Uncomment to initialize using linear profile instead of restart file
%     SHOULD REDO WITH ELLIPTIC SOLVE TO FIND INITIAL CONDITION FOR PHI
    z(nSpecies+1:nSpecies+1:end) = 0- + xCenter/xCenter(end)*numericalInfo.appliedVoltage;
    z(nSpecies+1) = 0;
    z(end) = numericalInfo.appliedVoltage;
end

% Output and plotting frequency. Will plot or produce output EVERY time
% step if the specified plotting/output time interval is less than the
% specified time step.
outputInt = max(floor(dtOut/dt), 1);
doOutput = numericalInfo.doOutput;
plotInt = max(floor(dtPlot/dt), 1);
doPlot = numericalInfo.doPlot;

%%%%% VOLTAGE SWEEP FUNCTION. TAKE ANOTHER LOOK AT THIS WITH HIGH PRECISION
%%%%% VARIABLES, BEFORE CODE IS PUBLISHED. NEED TO WRITE A BETTER
%%%%% DESCRIPTION OF THIS AND ENSURE THAT INPUT FORMAT IS CLEAR.
if numericalInfo.voltageSweepOn
    S = load(numericalInfo.voltageFileName);
    fieldList = fieldnames(S);
    voltageVals = getfield(S, fieldList{1});
    voltageRate = numericalInfo.voltageSweepRate;
    timeVals = tStart:(1/voltageRate):(length(voltageVals)-1)/voltageRate;

    if numericalInfo.highPrecision
        z(nSpecies+1:nSpecies+1:end) = mp(0 + voltageVals(1)*xCenter/xCenter(end));
        z(nSpecies+1) = mp(0);
        z(end) = mp(voltageVals(1));
    else
        z(nSpecies+1:nSpecies+1:end) = 0 + voltageVals(1)*xCenter/xCenter(end);
        z(nSpecies+1) = 0;
        z(end) = voltageVals(1);
    end
end

%output generation function here, to save initial condition
if doOutput
    produceOutput(xCenter, uniqueSpecies, z, t, folderName, numericalInfo.highPrecision);
end

% Plotting function here, to show initial condition
if doPlot
%         junctionThickness = numericalInfo.interfaces(4) - numericalInfo.interfaces(3);
%         centerPosition = 1/2*(numericalInfo.interfaces(4) + numericalInfo.interfaces(3));
        junctionThickness = 100e-6;
        centerPosition = 1330e-6;
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 23);
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 8);
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        %set(gca, 'YScale', 'log')
        %ylim([1e-16 2])
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'potential', 24);
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'Efield', 10);
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 25)
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 12)
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'waterEq', 26)
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'waterEq', 14)
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
        
%         producePlot(z, constants, xCenter, xFace, dxC, dxF, [2 5], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 7);
%         drawnow
end

% Initialization of checkpoint variables and counter. Every
% <checkPointInterval> successful solves, a snapshot of the current z, z at
% the previous timestep, and the previous dt are saved. If the adaptive
% time refinement requires refinement below <dtLowerThreshhold>, then the
% code jumps back to the checkpoint and attempts to solve from that point,
% using smaller dt values.
solveCounter = 0;
checkPointInterval = 16;
dtLowerThreshhold = 1e-16;
zCheckPoint = z;
zPreviousTimestepCheckpoint = z;
dtPreviousTimestepCheckPoint = -1;

% Initialization of tracker for effective time step in previous solution
% attempt. 
dtEffOld = dt;

% Initialization of dummy variables for first step - eventually used to
% store old profiles and time step, for second order temporal scheme.
dtPreviousTimestep = -1;
zPreviousTimestep = z;

% Create row and column index vectors for the sparse matrix used in the
% linear solve. This only needs to be done once, so we do it here outside
% the time stepping loop.
% [rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, COO2CSCIndices, A] = ...
%     makeSparseMatrixIndexVectors_fourthOrder(length(uniqueSpecies), length(xCenter), 10);
[rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, COO2CSCIndices, A] = ...
    makeSparseMatrixIndexVectors(length(uniqueSpecies), length(xCenter));

% Start timing, for completion estimate calculation
tic;

hasCheckpointedOnce = 0;

% Time-stepping loop
for ii = (ceil(t/dt)+1):totalSteps
    % Initialize convergence check, effective dt. Set target time.
    notConverged = 1;
    dtEff = dt;
    targetTime = t+dt;

    % Do time step with size dtEff if prior attempt did not converge or if
    % we are not yet at the target time.
    while (notConverged) || (t ~= targetTime)

        tEndEff = t + dtEff;

        % Adjust boundary condition for electric potential, sweep along
        % user-specified profile.
        if numericalInfo.voltageSweepOn
            if tEndEff <= timeVals(end)
                z(end) = interp1(timeVals, voltageVals, tEndEff);
            else
                z(end) = voltageVals(end);
            end
        end

        % PERFORM TIME STEP
        [zOut, tOut] = doTimeStepThermEq_highOrder(dxC, dxF, nSpecies, z, zPreviousTimestep, t, dtEff, dtPreviousTimestep, constants, rowInd, colInd, nnz, nEq, ...
            valsArrayBlockSizePerSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], numericalInfo.highPrecision, ...
            COO2CSCIndices, A);

        % Check for convergence: Are there any negative values for species
        % concentrations? Are there any NAN values or imaginare values
        % anywhere?
        negVals = 0;
        for jj = 1:nSpecies
            negVals = negVals || any(zOut(jj:(nSpecies+1):end)<-eps);
        end
        nanVals = any(isnan(zOut));
        imagVals = any(~isreal(zOut));

        if negVals ||nanVals || imagVals
            % Simulation is not converged, time refinement is required.
            notConverged = 1;
            % If larger than lower threshhold, halve timestep and try again
            if dtEff > dtLowerThreshhold
                % If refinement is required, ensure that <dtEff> is never
                % more than 2x larger than most recent successful dtEff
                % value (<dtEffOld>)
                dtEff = min([dtEff/2, 2*dtEffOld]);
                disp(['Smallest time step: ' num2str(dtEff)]);
            else
                % Revert to checkpoint and try again using <dtEff>/2
                if hasCheckpointedOnce
                    return
                end
                z = zCheckPoint;
                zPreviousTimestep = zPreviousTimestepCheckpoint;
                dtPreviousTimestep = dtPreviousTimestepCheckPoint;
                dtEff = min([dtEff/2, 2*dtEffOld]);
                disp(['Smallest time step (CHECKPOINT REQUIRED): ' num2str(dtEff)]);
                
                hasCheckpointedOnce = 1;
            end
        else
            % Simulation is converged, save successful dtEff in <dtEffOld>,
            % save previous time step data for second-order temporal
            % scheme.
            notConverged = 0;
            solveCounter = solveCounter + 1;
            zPreviousTimestep = z;
            z = zOut;
            t = tOut;
            dtEffOld = dtEff;
            dtPreviousTimestep = dtEff;

            % Increase the time-step if possible, but never more than 2x
            % the most recent successful dtEff.
            dtEff = min([targetTime - t, 2*dtEffOld]);

            % Save checkpoint data if necessary.
            if mod(solveCounter, checkPointInterval) == 0
                zCheckPoint = z;
                zPreviousTimestepCheckpoint = zPreviousTimestep;
                dtPreviousTimestepCheckPoint = dtPreviousTimestep;
            end
            
            disp(['Reached time: ' num2str(t)])
        end
    end

    if (mod(ii, outputInt) == 0) && (doOutput)
        % output generation function here
        produceOutput(xCenter, uniqueSpecies, z, t, folderName, numericalInfo.highPrecision);
        
        if ~ doPlot
            % Completion time estimate
            elapsedTimeSeconds = toc;
            completedFraction = ii/totalSteps;
            stepsBypassedDueToRestart = ceil(tStart/dt);
            estimatedTimeRemainingSeconds = (1-completedFraction) * elapsedTimeSeconds/((ii - stepsBypassedDueToRestart)/totalSteps);

            disp(['Completed step ' num2str(ii) ' of ' num2str(totalSteps)...
                ' (' num2str(100*completedFraction) '%) t = ' num2str(t) ...
                '. Est. remaining: ' num2str(estimatedTimeRemainingSeconds/60/60) ' hr.'])
        end
    end

    if (mod(ii, plotInt) == 0) && doPlot
%         junctionThickness = numericalInfo.interfaces(4) - numericalInfo.interfaces(3);
%         centerPosition = 1/2*(numericalInfo.interfaces(4) + numericalInfo.interfaces(3));
        junctionThickness = 100e-6;
        centerPosition = 1330e-6;
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 23);
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 8);
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        %set(gca, 'YScale', 'log')
        %ylim([1e-16 2])
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'potential', 24);
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'Efield', 10);
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 25)
        drawnow
%         
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 2 3 4 5 6 7], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 12)
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'waterEq', 26)
        drawnow
        
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'waterEq', 14)
        xlim([centerPosition - junctionThickness, centerPosition + junctionThickness])
        drawnow
        
%         producePlot(z, constants, xCenter, xFace, dxC, dxF, [2 5], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 7);
%         drawnow

        % Completion time estimate
        elapsedTimeSeconds = toc;
        completedFraction = ii/totalSteps;
        stepsBypassedDueToRestart = ceil(tStart/dt);
        estimatedTimeRemainingSeconds = (1-completedFraction) * elapsedTimeSeconds/((ii - stepsBypassedDueToRestart)/totalSteps);

        disp(['Completed step ' num2str(ii) ' of ' num2str(totalSteps)...
            ' (' num2str(100*completedFraction) '%) t = ' num2str(t) ...
            '. Est. remaining: ' num2str(estimatedTimeRemainingSeconds/60/60) ' hr.'])
    end
end
end
%------------- END OF CODE --------------
