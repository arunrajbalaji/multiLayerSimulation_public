% Simulation to determine chemical equilibrium for a certain layer. Uses
% input that was generated in the generateChemEqInputFile() method. Runs
% a simulation out until steady state, in order to determine initial
% condition for the target simulation.
function equilibratedValues = runChemEqSimulation(folderName, fileName)
% Input folder and file name
fullName = [folderName, fileName];

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
[xCenter, xFace, dxC, dxF] = genOverallMeshConstantEDLRes(numericalInfo.interfaces,...
    [layerInfo(:).dxMax], [layerInfo(:).dxMin],...
    [layerInfo(:).LID], [layerInfo(:).gridSymmetry], [0, 0, 0, 0, 0]);
% [xCenter, xFace, dxC, dxF] = genOverallMeshConstantEDLRes(numericalInfo.interfaces,...
%     [layerInfo(:).dxMax], [layerInfo(:).dxMin],...
%     [layerInfo(:).LID], [layerInfo(:).gridSymmetry], [1e-9]);

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
        z((end-nSpecies):1:end-1) = constants.initVal(end,:);
    end
%    z(nSpecies+1:nSpecies+1:end) = interp1(xCenter_old, rawData(nSpecies+2,:), xCenter);
%     Uncomment to initialize using linear profile instead of restart file
%     SHOULD REDO WITH ELLIPTIC SOLVE TO FIND INITIAL CONDITION FOR PHI
    z(nSpecies+1:nSpecies+1:end) = 0 + xCenter/xCenter(end)*numericalInfo.appliedVoltage;
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
    producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 1);
%     ylim([0 1.25])
    drawnow
    %     xlim([4e-2-100e-6 4e-2+230e-6+100e-6])
    %     drawnow
    producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'potential', 2);
    %     xlim([4e-2-100e-6 4e-2+230e-6+100e-6])
    %         drawnow
    %     producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'flux', 3)
    producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 4)
    producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'charge', 5);
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
[rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, COO2CSCIndices, A] = ...
    makeSparseMatrixIndexVectors(length(uniqueSpecies), length(xCenter));

% Start timing, for completion estimate calculation
tic;

% Compute original charge, for comparison later and determination of
% whether or not failure has occurred.
% chargeProfile = constants.backCharge * constants.nA * constants.e * constants.litersPerCubicMeter;
% for speciesIndex = 1:nSpecies
%     chargeProfile = chargeProfile + constants.vale(1,speciesIndex) ...
%         * z(speciesIndex:(nSpecies+1):end)  * constants.nA * constants.e * constants.litersPerCubicMeter;
% end
% originalTotalCharge = sum(dxC .* chargeProfile);

% original value, for comparison to compute delta
zPrevious = z;

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
        [zOut, tOut] = doTimeStep(dxC, dxF, nSpecies, z, zPreviousTimestep, t, dtEff, dtPreviousTimestep, constants, rowInd, colInd, nnz, nEq, ...
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

%         chargeProfile = constants.backCharge * constants.nA * constants.e * constants.litersPerCubicMeter;
%         for speciesIndex = 1:nSpecies
%             chargeProfile = chargeProfile + constants.vale(1,speciesIndex) ...
%                 * zOut(speciesIndex:(nSpecies+1):end)  * constants.nA * constants.e * constants.litersPerCubicMeter;
%         end
%         totalCharge = sum(dxC .* chargeProfile);
% 
%         if abs(totalCharge - originalTotalCharge) > 100*eps
%             chargeConservationFailure = 1;
%             abs(totalCharge - originalTotalCharge)
%         else
%             chargeConservationFailure = 0;
%         end

        chargeConservationFailure = 0;

        if negVals ||nanVals || imagVals || chargeConservationFailure
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
                z = zCheckPoint;
                zPreviousTimestep = zPreviousTimestepCheckpoint;
                dtPreviousTimestep = dtPreviousTimestepCheckPoint;
                dtEff = min([dtEff/2, 2*dtEffOld]);
                disp(['Smallest time step (CHECKPOINT REQUIRED): ' num2str(dtEff)]);
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

%             disp(['Reached time: ' num2str(t)])
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
        % Plotting function here
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'concentration', 1);
%         ylim([0 2.5])
%         xlim([0 200005e-9])
        drawnow
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'potential', 2);
%         xlim([0 200005e-9])
        drawnow
    %     producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'flux', 3)
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [1 3 4 5 6], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'current', 4)
        producePlot(z, constants, xCenter, xFace, dxC, dxF, [], uniqueSpecies, [numericalInfo.leftElectrodeBC, numericalInfo.rightElectrodeBC], 'charge', 5);

        % Completion time estimate
        elapsedTimeSeconds = toc;
        completedFraction = ii/totalSteps;
        stepsBypassedDueToRestart = ceil(tStart/dt);
        estimatedTimeRemainingSeconds = (1-completedFraction) * elapsedTimeSeconds/((ii - stepsBypassedDueToRestart)/totalSteps);

        disp(['Completed step ' num2str(ii) ' of ' num2str(totalSteps)...
            ' (' num2str(100*completedFraction) '%) t = ' num2str(t) ...
            '. Est. remaining: ' num2str(estimatedTimeRemainingSeconds/60/60) ' hr.'])
    end
    
   maxVal = 0;
   for sIndex = 1:nSpecies
        maxVal = max(max(abs(z(sIndex:(nSpecies+1):end))-zPrevious(sIndex:(nSpecies+1):end)), maxVal);
   end
   if mod(ii,plotInt) == 0
        disp(['Chemical Equilibrium Solve Max Error: ' num2str(maxVal)])
   end
   zPrevious = z;
end

 equilibratedValues = z(2*(nSpecies+1)+1:3*(nSpecies+1)-1);
end