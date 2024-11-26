function [z, t] = doTimeStepThermEq_highOrder(dxC, dxF, nSpecies, z, zPreviousTimeStep, t, dt, ...
    dtPreviousTimeStep, constants, rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, ...
    isElectrodeBC, highPrecision, COO2CSCIndices, A)
% doTimeStep: Perform the time step using first order implicit in time for
% the first time step, second order implicit in time for all subsequent
% time steps, and second-order central in space. Linear operator formation
% must be performed in the same order as specified by rowInd and coldInd,
% which are formed in makeSparseMatricIndexVectors.m. Constants are read in
% using the struct formed buy genPhysConstArrays.m and parseInputFile.m.
%
%   [z, t] = doTimeStep(dxC, dxF, nSpecies, z, zPreviousTimeStep, t, dt, ...
%    dtPreviousTimeStep, constants, rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, ...
%    highPrecision)
%
% Inputs:
%       dxC             - Mesh size at cell centers.
%       dxF             - Mesh size at cell faces.
%       nSpecies        - Number of unique species
%       z               - Unknowns vector before step is performed. See
%                           initDomain.m for details of unknowns storage.
%       zPreviousTimeStep 
%                        - Unknowns vector at end of previous time step.
%       t               - time before step is performed.
%       dt              - time step to advance with
%       dtPreviousTimeStep
%                       - Time step used between z and zPreviousTimeStep.
%                         Set equal to -1 if you would like to use
%                         first-order scheme in time (i.e. for very first
%                         times step, when no previous steps are available)
%       constants       - physical constants, in struct as per output of
%                         genPhysConstArrays.m and parseInputFile.m
%       rowInd          - Sparse system row indices. See
%                         makeSparseMatrixIndexVectors.m
%       colInd          - Sparse system row indices. See
%                         makeSparseMatrixIndexVectors.m
%       nnz             - Number of non-zeros in the matrix
%       nEq             - Number of unknowns to solve for. Remember that
%                         first and last cell are ghost cells, not solved.
%       valsArrayBlockSizePerSpecies
%                       - Number of non-zero entries per species. Useful
%                         for formation of vals vector.
%       isElectrodeBC   - 1 for electrode boundary condition, 0 for
%                         reservoir (Dirichlet) boundary condition. Vector
%                         of size two, representing left and right side
%                         respectively.
%       highPrecision   - FLAG, 1 for advanpix high-precision mode,
%
% Outputs:
%       z               - Unknows vector after time step
%       t               - Time after time step
%
% Other m-files required: none
% MAT-files required: none
%
% See also: parseInputFile.m, genPhysConstArrays.m,
%           initDomain.m, makeSparseMatrixIndexVectors.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2021
%------------- BEGIN CODE --------------
% Unpack constants
poro = constants.poro;
tort = constants.tort;
perm = constants.perm;
diff = constants.diff;
acti = constants.acti;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
nA = constants.nA;
litersPerCubicMeter = constants.litersPerCubicMeter;
backCharge = constants.backCharge;
wienBeta = constants.wienBeta;
leftElectrodeBC = isElectrodeBC(1);
rightElectrodeBC = isElectrodeBC(2);

% wienActivityExponent = [-1/2 0 -0 0 0 -1/2 0];
%%%%% USE THIS ONE
%wienActivityExponent = [-1/2 0 0 0 -1/2 -1/2 -1/2];

% wienActivityExponent = [0 0 0 0 0 0 0];

wienActivityExponent = constants.wienActivityExponent;

rxnRateCoefficients = constants.rxnRateCoefficients;
rxnInputOutputIndices = constants.rxnInputOutputIndices;
rxnWienCoeff = constants.rxnWienCoeff;

deltaStern = constants.deltaStern;
leftExchCurrent = constants.leftIonExchangeCurrent;
rightExchCurrent = constants.rightIonExchangeCurrent;

faradaicCoeff = deltaStern*e/kb/T;

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

if highPrecision
    vals = mp(zeros(nnz,1));
    b = mp(zeros(nEq,1));
else
    vals = zeros(nnz,1);
    b = zeros(nEq,1);
end

t = t + dt;
N = length(dxC);

% Loop over iterations
nIter = 3;
zOriginal = z;

% ASSUMING SINGLE VALUE OF WIEN COEFFICIENT, GIVEN BY LARGEST VALUE
wienAlphaBeta = max(rxnWienCoeff)*wienBeta;

% Ad hoc parameters for second wien higher order implementation
% Applied for all points except for outermost points, so that boundary
% condition handling is not required. Should be fine since electric field
% is weak in these regions.
dPhiLeftLeftAll = zeros(nSpecies*(N-6), 1);
dPhiRightRightAll = zeros(nSpecies*(N-6), 1);
rowIndLeftAll = zeros(nSpecies*(N-6), 1);
rowIndRightAll = zeros(nSpecies*(N-6), 1);
colIndLeftAll = zeros(nSpecies*(N-6), 1);
colIndRightAll = zeros(nSpecies*(N-6), 1);

for ii = 1:nIter
    %% Determine LHS matrix A values and RHS vector b values
    % Choose between first and second order scheme in time
    if dtPreviousTimeStep < 0
        timeTerm = poro(2:end-1)/dt;
    else
        timeGamma = 1/(dtPreviousTimeStep + dtPreviousTimeStep^2/dt);
        timeBeta = -timeGamma*(dt + dtPreviousTimeStep)^2/dt^2;
        timeAlpha = -timeGamma - timeBeta;
        timeTerm = poro(2:end-1) * timeAlpha;
    end
    
    % Compute summation term for steric effect flux.
    reshapedZ = reshape(z, nSpecies+1, length(dxC));
    stericTermSpeciesSum = transpose(sum(reshapedZ(1:nSpecies, :).*stericOnOffVec, 1));
    
    % dPhidX compted at cell centers, size = (N-2). Assumes Dirichlet,
    % useful for reaction term (second Wien effect) - exponential form
    dPhidXStar = (z(3*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1)))./(dxF(1:end-1) + dxF(2:end));
    
    % Outermost species loop - used to assemble matrix entries in vals,
    % vectorized for each species.
    for jj = 1:nSpecies
        % outer derivative, length = (N-2) associated with cell centers
        % that are actively computed
        outerDerivativeRightFace = 1./dxC(2:end-1).*(poro(2:end-1)+poro(3:end))/2.*(diff(2:end-1,jj)+diff(3:end,jj))/2./((tort(2:end-1)+tort(3:end))/2);
        outerDerivativeLeftFace = 1./dxC(2:end-1).*(poro(1:end-2)+poro(2:end-1))/2.*(diff(1:end-2,jj)+diff(2:end-1,jj))/2./((tort(1:end-2)+tort(2:end-1))/2);
        
        % Diffusion flux, length = (N-1), computed at faces (left and
        % right) of all cells that are actively computed. Considered
        % 'incomplete' because of need to multiply by activity coeff again.
        matrixDiffusionFluxIncomplete = 1./dxF.*2./(acti(1:end-1, jj) + acti(2:end, jj));
        if leftElectrodeBC
            matrixDiffusionFluxIncomplete(1) = 0;
        end
        if rightElectrodeBC
            matrixDiffusionFluxIncomplete(end) = 0;
        end
        
        % Electromigration flux, length = (N-1), computed at faces (left
        % and right) of all cells that are actively computed.
        matrixElectromigrationFlux = 1/2*(vale(1,jj)*e/kb/T)*(z(2*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-(nSpecies+1)))./dxF;
        if leftElectrodeBC
            matrixElectromigrationFlux(1) = 0;
        end
        if rightElectrodeBC
            matrixElectromigrationFlux(end) = 0;
        end
        
        % Steric flux, length = (N-1), computed at faces (left and right)
        % of all cells that are actively computed
        matrixStericFlux = stericOnOffVec(jj) * 1/2 * log((1 - stericACubed * stericTermSpeciesSum(2:end))./(1 - stericACubed * stericTermSpeciesSum(1:end-1)))./dxF;
        if leftElectrodeBC
            matrixStericFlux(1) = 0;
        end
        if rightElectrodeBC
            matrixStericFlux(end) = 0;
        end
        
        % Steric effect cross term, length = (N-1), computed at faces (left
        % and right) of all cells that are actively computed. Incomplete
        % because this is missing the a^3/(1 - a^3 * SUM), which is
        % multiplied in later.
        matrixStericCrossTermIncomplete = stericOnOffVec(jj).*(z(jj:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+jj:(nSpecies+1):end))/2 ./dxF;
        if leftElectrodeBC
            matrixStericCrossTermIncomplete(1) = 0;
        end
        if rightElectrodeBC
            matrixStericCrossTermIncomplete(end) = 0;
        end
        
        % Coefficient for dPhi terms in electromigration, length = (N-1),
        % computed at faces of all cells that are actively computed.
        matrixElectromigrationPotentialTerm = vale(1,jj)*e/kb/T*(z(jj:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+jj:(nSpecies+1):end))/2./dxF;
        if leftElectrodeBC
            matrixElectromigrationPotentialTerm(1) = 0;
        end
        if rightElectrodeBC
            matrixElectromigrationPotentialTerm(end) = 0;
        end
        
        % Implicit terms for second Wien effect (ignoring boundary cells)
        wienFieldStrength = wienAlphaBeta .* abs(z(3*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1)))./(dxF(2:end) + dxF(1:end-1));
        matrixWienTermDeltaC = wienActivityExponent(jj) * (1./dxF(2:end-1)) .* (wienFieldStrength(2:end) - wienFieldStrength(1:end-1));
        
        cInterp = 1/2 * (z((nSpecies+1)+jj:(nSpecies+1):end) + z(jj:(nSpecies+1):end-(nSpecies+1)));
        dxC_double = dxF(1:end-1) + dxF(2:end);
        matrixWienTermDeltaPhi = wienActivityExponent(jj) * cInterp ./ dxF;
        
        % dc center
        dcCenterSpeciesTransport = timeTerm ...
                    + outerDerivativeRightFace.*(matrixDiffusionFluxIncomplete(2:end).*acti(2:end-1, jj) - matrixElectromigrationFlux(2:end) + matrixStericFlux(2:end)) ...
                    + outerDerivativeLeftFace.*(matrixDiffusionFluxIncomplete(1:end-1).*acti(2:end-1, jj) + matrixElectromigrationFlux(1:end-1) - matrixStericFlux(1:end-1)) ...
                    + 1/2*[0; 0; -outerDerivativeRightFace(3:end-2).*matrixWienTermDeltaC(3:end-1) + outerDerivativeLeftFace(3:end-2).*matrixWienTermDeltaC(2:end-2); 0; 0];
        
        % dc left
        dcLeftSpeciesTransport = outerDerivativeLeftFace(2:end).*(-matrixDiffusionFluxIncomplete(2:end-1).*acti(2:end-2,jj) + matrixElectromigrationFlux(2:end-1)  - matrixStericFlux(2:end-1)) ...
            + 1/2 * [0; outerDerivativeLeftFace(3:end-2) .* matrixWienTermDeltaC(2:end-2); 0; 0];
        
        % dc right
        dcRightSpeciesTransport = outerDerivativeRightFace(1:end-1).*(-matrixDiffusionFluxIncomplete(2:end-1).*acti(3:end-1,jj) - matrixElectromigrationFlux(2:end-1)  + matrixStericFlux(2:end-1)) ...
            - 1/2 * [0; 0; outerDerivativeRightFace(3:end-2).*matrixWienTermDeltaC(3:end-1); 0];
        
        % dc cross terms center for steric term
        dcCenterStericCrossTerms = repmat((stericACubed./(1 - stericACubed * stericTermSpeciesSum(2:end-1))).*(outerDerivativeRightFace.*matrixStericCrossTermIncomplete(2:end) + outerDerivativeLeftFace.*matrixStericCrossTermIncomplete(1:end-1)), [nSpecies, 1]);
        
        % dc cross terms left for steric term
        dcLeftStericCrossTerms = repmat(-outerDerivativeLeftFace(2:end) .* matrixStericCrossTermIncomplete(2:end-1) .* (stericACubed./(1 - stericACubed * stericTermSpeciesSum(2:end-2))), [nSpecies, 1]);
        
        % dc cross terms right for steric term
        dcRightStericCrossTerms = repmat(-outerDerivativeRightFace(1:end-1) .* matrixStericCrossTermIncomplete(2:end-1) .* (stericACubed./(1 - stericACubed * stericTermSpeciesSum(3:end-1))), [nSpecies, 1]);
        
        % dPhi center
        dphiCenterSpeciesTransport = outerDerivativeRightFace.*matrixElectromigrationPotentialTerm(2:end) + outerDerivativeLeftFace.*matrixElectromigrationPotentialTerm(1:end-1) ...
            + [0; 0; sign(dPhidXStar(4:end-1)) .* (outerDerivativeRightFace(3:end-2) .* matrixWienTermDeltaPhi(4:end-2) .* wienAlphaBeta ./ dxC_double(4:end-1)) ...
            - sign(dPhidXStar(2:end-3)) .* (outerDerivativeLeftFace(3:end-2) .* matrixWienTermDeltaPhi(3:end-3) .* wienAlphaBeta ./ dxC_double(2:end-3)); 0; 0];
        
        % dPhi left
        dphiLeftSpeciesTransport = - outerDerivativeLeftFace(2:end) .* matrixElectromigrationPotentialTerm(2:end-1) ...
            + [0; -sign(dPhidXStar(3:end-2)) .* (outerDerivativeLeftFace(3:end-2) .* matrixWienTermDeltaPhi(3:end-3) .* wienAlphaBeta ./ dxC_double(3:end-2)) ...
            - sign(dPhidXStar(3:end-2)) .* (outerDerivativeRightFace(3:end-2) .* matrixWienTermDeltaPhi(4:end-2) .* wienAlphaBeta ./ dxC_double(3:end-2)); 0; 0];
        
        % dPhi right
        dphiRightSpeciesTransport = -outerDerivativeRightFace(1:end-1) .* matrixElectromigrationPotentialTerm(2:end-1) ...
            + [0; 0; sign(dPhidXStar(3:end-2)) .* (outerDerivativeLeftFace(3:end-2) .* matrixWienTermDeltaPhi(3:end-3) .* wienAlphaBeta ./ dxC_double(3:end-2)) ...
            + sign(dPhidXStar(3:end-2)) .* (outerDerivativeRightFace(3:end-2) .* matrixWienTermDeltaPhi(4:end-2) .* wienAlphaBeta ./ dxC_double(3:end-2)); 0];
        
        % Wider stencil for wien effect activity term
        % For now, just adding to the end of the entire vector.
        % Can be repackaged for speed if necessary
        dPhiRightRightAll(jj:nSpecies:end) = -sign(dPhidXStar(4:end-1)) .* (outerDerivativeRightFace(3:end-2) .* matrixWienTermDeltaPhi(4:end-2) .* wienAlphaBeta ./ dxC_double(4:end-1));
        dPhiLeftLeftAll(jj:nSpecies:end) = + sign(dPhidXStar(2:end-3)) .* (outerDerivativeLeftFace(3:end-2) .* matrixWienTermDeltaPhi(3:end-3) .* wienAlphaBeta ./ dxC_double(2:end-3));
        rowIndRightAll(jj:nSpecies:end) = (2*(nSpecies+1)+jj:(nSpecies+1):(N-4)*(nSpecies+1));
        rowIndLeftAll = rowIndRightAll;
        colIndRightAll(jj:nSpecies:end) = (5*(nSpecies+1):(nSpecies+1):(N-2)*(nSpecies+1));
        colIndLeftAll(jj:nSpecies:end) = ((nSpecies+1):(nSpecies+1):(N-6)*(nSpecies+1));

        % Allocate space for accumulation of coefficients due to reaction
        % terms and transport combined, for cross-species concentration
        % coupling and for elecric potential (Due to second Wien effect)
        dcCenterAllCrossTerms = zeros(size(dcCenterStericCrossTerms));
        dphiLeftTotal = zeros(size(dphiLeftSpeciesTransport));
        dphiRightTotal = zeros(size(dphiRightSpeciesTransport));
        RHSReactionTerms = zeros(length(dxC)-2, 1);
        
        % Loop over reactions to include reaction terms for species jj.
        for ll = 1:size(rxnInputOutputIndices, 1)
            
            % Set all reaction term accumulation vectors to zero initially.
            % Otherwise will have problems - variable scope is not
            % restricted to single iteration of each reaction loop.
            dcCenterReactionCrossTermsReactants = zeros(size(dcCenterStericCrossTerms));
            dcCenterReactionCrossTermsProductsPos3 = zeros(size(dcCenterStericCrossTerms));
            dcCenterReactionCrossTermsProductsPos4 = zeros(size(dcCenterStericCrossTerms));
            dphiLeftReactants = zeros(size(dphiLeftSpeciesTransport));
            dphiLeftProductsPos3 = zeros(size(dphiLeftSpeciesTransport));
            dphiLeftProductsPos4 = zeros(size(dphiLeftSpeciesTransport));
            dphiRightReactants = zeros(size(dphiRightSpeciesTransport));
            dphiRightProductsPos3 = zeros(size(dphiLeftSpeciesTransport));
            dphiRightProductsPos4 = zeros(size(dphiLeftSpeciesTransport));
            
            % CONDITIONAL: (if) Is it in reactant position 1
            if jj == rxnInputOutputIndices(ll, 1)
                % SUBCONDITIONAL: (if)Is it also in reactant position 2
                if jj == rxnInputOutputIndices(ll, 2)
                    dcRxnEq1Spec1 = 4*poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1));
                    
                    % Add reactant terms in Eq 1 (same as Eq 2) for Spec
                    % 1 (same as Spec 2) 
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq1Spec1;

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms - 2 * poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)).^2;
                    
                % SUBCONDITIONAL: (elseif) Is there a different species in reactant
                % position 2
                elseif rxnInputOutputIndices(ll, 2) ~= 0
                    % Add reactant terms in Equation 1 alone, Equation 2
                    % terms will be added when outer loop reaches the
                    % species index for that species (***)
                    dcRxnEq1Spec1 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1));
                    dcRxnEq1Spec2 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                    
                    % Add reactant terms in Eq 1 for Spec 1
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq1Spec1;
                    
                    % Add reactant terms in Eq 1 for Spec 2
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq1Spec2;
                
                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms - poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));

                % SUBCONDITIONAL: (else) Is reactant position 2 empty?
                else
                    % Add reactant terms in Equation 1 alone, for Spec 1
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
%                     disp(['Max Wien exponential:' num2str(max(wienExponential))])
                    dcRxnEq1Spec1 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq1Spec1;
                    
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftReactants = -dphiTermForReactions(2:end);
                    dphiRightReactants = dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms - poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1));
                end
            % CONDITIONAL: (elseif) Is it in reactant position 2
            elseif jj == rxnInputOutputIndices(ll, 2)
                % SUBCONDITIONAL: (if) Is there a different species in reactant
                % position 1. (Cannot be same species, due to previous
                % conditional statement)
                if rxnInputOutputIndices(ll, 1) ~= 0
                    % Add reactant terms in Equation 2 alone. This also
                    % completes the equation above (***)
                    dcRxnEq2Spec1 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1));
                    dcRxnEq2Spec2 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                    
                    % Add reactant terms in Eq 2 for Spec 1
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq2Spec1;
                    
                    % Add reactant terms in Eq 2 for Spec 2
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq2Spec2;

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms - poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));
                    
                % SUBCONDITIONAL: (else) Is reactant position 1 empty
                else
                    % Add reactant terms to Equation 2 alone, for Spec 2
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    dcRxnEq2Spec2 = poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    dcCenterReactionCrossTermsReactants((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq2Spec2;
                    
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftReactants = -dphiTermForReactions(2:end);
                    dphiRightReactants = dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms - poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));
                    
                end
            end
            
            % CONDITIONAL: (if) Is it in product position 3
            if jj == rxnInputOutputIndices(ll, 3)
                % SUBCONDITIONAL: (if) there are reactants in position 1
                % and position 2
                if rxnInputOutputIndices(ll,1)~=0 && rxnInputOutputIndices(ll,2)~=0
                    if rxnInputOutputIndices(ll,1) == rxnInputOutputIndices(ll,2)
                        dcRxnEq3Spec1 = -2 * poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                    
                        % Add term in Eq 3 for species 1, same as species 2
                        dcCenterReactionCrossTermsProductsPos3((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq3Spec1;
                        
                        % Add reaction term to right hand side
                        RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)).^2;
                    else
                        dcRxnEq3Spec1 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1));
                        dcRxnEq3Spec2 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                    
                        % Add term in Eq 3 for species 1
                        dcCenterReactionCrossTermsProductsPos3((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq3Spec1;
                        
                        % Add term in Eq 3 for species 2
                        dcCenterReactionCrossTermsProductsPos3((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq3Spec2;
                        
                        % Add reaction term to right hand side
                        RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));
                    end

                % SUBCONDITIONAL: (elseif) there are reactants only in
                % position 1
                elseif rxnInputOutputIndices(ll,1) ~= 0
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    dcRxnEq3Spec1 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    % Add term in Eq 3 for species 1
                    dcCenterReactionCrossTermsProductsPos3((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq3Spec1;
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftProductsPos3 = dphiTermForReactions(2:end);
                    dphiRightProductsPos3 = -dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1));

                % SUBCONDITIONAL: (elseif) there are reactants only in
                % position 2
                elseif rxnInputOutputIndices(ll,2) ~= 0
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    dcRxnEq3Spec2 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    % Add term in Eq 3 for species 2
                    dcCenterReactionCrossTermsProductsPos3((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq3Spec2;
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftProductsPos3 = dphiTermForReactions(2:end);
                    dphiRightProductsPos3 = -dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));

                % SUBCONDITION: (else) there are no reactants, only a
                % constant source term
                else
                    % Add reaction term to right hand side
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential;
                end
            end
            
            % CONDITIONAL: (if) Is it in product position 4
            if jj == rxnInputOutputIndices(ll, 4)
                % SUBCONDITIONAL: (if) there are reactants in position 1
                % and position 2
                if rxnInputOutputIndices(ll,1)~=0 && rxnInputOutputIndices(ll,2)~=0
                    if rxnInputOutputIndices(ll,1) == rxnInputOutputIndices(ll,2)
                        dcRxnEq4Spec1 = -2*poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                    
                        % Add term in Eq 4 for species 1, same as species 2
                        dcCenterReactionCrossTermsProductsPos4((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq4Spec1;
                        
                        % Add reaction term to right hand side
                        RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)).^2;
                    else
                        dcRxnEq4Spec1 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1));
                        dcRxnEq4Spec2 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1));
                        
                        % Add term in Eq 3 for species 1
                        dcCenterReactionCrossTermsProductsPos4((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq4Spec1;
                        
                        % Add term in Eq 3 for species 2
                        dcCenterReactionCrossTermsProductsPos4((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq4Spec2;
    
                        % Add reaction term to right hand side
                        RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1)) .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));
                    end
                    
                % SUBCONDITIONAL: (elseif) there are reactants only in
                % position 1
                elseif rxnInputOutputIndices(ll,1) ~= 0
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    dcRxnEq4Spec1 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    % Add term in Eq 4 for species 1
                    dcCenterReactionCrossTermsProductsPos4((rxnInputOutputIndices(ll, 1)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 1))*(length(dxC)-2)) = dcRxnEq4Spec1;
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,1):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftProductsPos4 = dphiTermForReactions(2:end);
                    dphiRightProductsPos4 = -dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 1):(nSpecies+1):end-(nSpecies+1));
                    
                % SUBCONDITIONAL: (elseif) there are reactants only in
                % position 2
                elseif rxnInputOutputIndices(ll,2) ~= 0
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    dcRxnEq4Spec2 = -poro(2:end-1).*rxnRateCoefficients(2:end-1, ll).*wienExponential;
                    
                    % Add term in Eq 4 for species 2
                    dcCenterReactionCrossTermsProductsPos4((rxnInputOutputIndices(ll, 2)-1)*(length(dxC)-2)+1:(rxnInputOutputIndices(ll, 2))*(length(dxC)-2)) = dcRxnEq4Spec2;
                    dphiTermForReactions = poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* rxnWienCoeff(ll) .* wienBeta .* z((nSpecies+1)+rxnInputOutputIndices(ll,2):(nSpecies+1):end-(nSpecies+1)) .* sign(dPhidXStar) ./(2*dxC(2:end-1));
                    dphiLeftProductsPos4 = dphiTermForReactions(2:end);
                    dphiRightProductsPos4 = -dphiTermForReactions(1:end-1);

                    % Add reaction term to right hand side
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential .* z((nSpecies+1)+rxnInputOutputIndices(ll, 2):(nSpecies+1):end-(nSpecies+1));
                    
                % SUBCONDITION: (else) there are no reactants, only a
                % constant source term
                else
                    % Add reaction term to right hand side
                    wienExponential = exp(rxnWienCoeff(ll)*wienBeta*abs(dPhidXStar));
                    RHSReactionTerms = RHSReactionTerms + poro(2:end-1) .* rxnRateCoefficients(2:end-1, ll) .* wienExponential;
                end
            end
            
            dcCenterAllCrossTerms = dcCenterAllCrossTerms + dcCenterReactionCrossTermsReactants ...
                + dcCenterReactionCrossTermsProductsPos3 + dcCenterReactionCrossTermsProductsPos4;
            dphiLeftTotal = dphiLeftTotal + dphiLeftReactants ...
                + dphiLeftProductsPos3 + dphiLeftProductsPos4;
            dphiRightTotal = dphiRightTotal + dphiRightReactants ...
                + dphiRightProductsPos3 + dphiRightProductsPos4;
        end
        
        dcCenterAllCrossTerms = dcCenterAllCrossTerms + dcCenterStericCrossTerms ...
            + [zeros((jj-1)*(length(dxC)-2),1); dcCenterSpeciesTransport; zeros((nSpecies-jj)*(length(dxC)-2),1)];
        dcLeftTotal = dcLeftStericCrossTerms ...
            + [zeros((jj-1)*(length(dxC)-3),1); dcLeftSpeciesTransport; zeros((nSpecies-jj)*(length(dxC)-3),1)];
        dcRightTotal = dcRightStericCrossTerms ...
            + [zeros((jj-1)*(length(dxC)-3),1); dcRightSpeciesTransport; zeros((nSpecies-jj)*(length(dxC)-3),1)];
        dphiLeftTotal = dphiLeftTotal + dphiLeftSpeciesTransport;
        dphiRightTotal = dphiRightTotal + dphiRightSpeciesTransport;
        
        % Accumulate all values into vals matrix
        vals((jj-1)*(valsArrayBlockSizePerSpecies)+1:jj*(valsArrayBlockSizePerSpecies)) = ...
            [dcCenterAllCrossTerms; ...
            dcLeftTotal; ...
            dcRightTotal; ...
            dphiCenterSpeciesTransport; ...
            dphiLeftTotal; ...
            dphiRightTotal];

        % Determine righthand side vector, b, for each species transport
        % equation, using terms computed above as able to.
        
        % RH Time Term
        if dtPreviousTimeStep < 0
            timeTermRHS = -poro(2:end-1).*(z((nSpecies+1)+jj:(nSpecies+1):end-(nSpecies+1))-zOriginal((nSpecies+1)+jj:(nSpecies+1):end-(nSpecies+1)))/dt;
        else
            timeTermRHS = -poro(2:end-1).*(timeAlpha*z((nSpecies+1)+jj:(nSpecies+1):end-(nSpecies+1)) ...
                + timeBeta*zOriginal((nSpecies+1)+jj:(nSpecies+1):end-(nSpecies+1)) ...
                + timeGamma*zPreviousTimeStep((nSpecies+1)+jj:(nSpecies+1):end-(nSpecies+1)));
        end

        % RHS Diffusion Term. Length = N-1, since computed on cell faces
        % and then differentiated onto cell centers when b is computed.
        RHSDiffusionFlux = matrixDiffusionFluxIncomplete .* (z((nSpecies+1)+jj:(nSpecies+1):end).*acti(2:end,jj) - z(jj:(nSpecies+1):end-(nSpecies+1)).*acti(1:end-1, jj));

        % RHS Electromigration Term. Length = N-1, since computed on cell
        % faces and then differentiated onto cell centers when b is
        % computed.
        RHSElectromigrationFlux = matrixElectromigrationFlux .* (z(jj:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+jj:(nSpecies+1):end));

        % RHS Steric Effect Term. Length N-1, since computed on cell faces
        % and then differentiated onto cell centers when b is computed.
        RHSStericFlux = -matrixStericFlux .* (z(jj:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+jj:(nSpecies+1):end));
        
        % Faradaic reaction terms. Only need to add values to first and
        % last cell as necessary for each species. Ion exchange current
        % denity of zero -> no Faradaic reactions.
        if vale(1,jj) ~= 0      % vale~=0, therwise might divide by zero
            exchangeCurrentConversion = 1/(vale(1,jj) * constants.e) * 1/constants.nA * 1/constants.litersPerCubicMeter;
            leftFaradaicFlux = -exchangeCurrentConversion*leftExchCurrent(jj)*2*sinh(faradaicCoeff*(z(2*(nSpecies+1)) - z(nSpecies+1))/dxF(1))/dxC(2);
            rightFaradaicFlux = exchangeCurrentConversion*rightExchCurrent(jj)*2*sinh(faradaicCoeff*(z(end) - z(end-(nSpecies+1)))/dxF(end))/dxC(end-1);
            RHSFaradaicReactions = [leftFaradaicFlux; zeros(length(dxC)-4, 1); rightFaradaicFlux];
        else
            RHSFaradaicReactions = zeros(length(dxC)-2, 1);
        end
        
        RHSGammaE = matrixWienTermDeltaC .* cInterp(2:end-1);
        
        b(jj:(nSpecies+1):end) = timeTermRHS ...
            + outerDerivativeRightFace.*(RHSDiffusionFlux(2:end) + RHSElectromigrationFlux(2:end) + RHSStericFlux(2:end)) ...
            - outerDerivativeLeftFace.*(RHSDiffusionFlux(1:end-1) + RHSElectromigrationFlux(1:end-1) + RHSStericFlux(1:end-1)) ...
            + RHSReactionTerms + RHSFaradaicReactions ...
            + [0; 0; (outerDerivativeRightFace(3:end-2) .* RHSGammaE(3:end-1) ...
            - outerDerivativeLeftFace(3:end-2) .* RHSGammaE(2:end-2)); 0; 0];
    end
    
    % Gauss's Equation: coupling to center, left, right phi terms and to
    % all local species terms

    % Permittivity * dphi/dx, length = (N-1), computed at faces (left and
    % right) of all cells that are actively computed.
    permDphiDx = (perm(2:end)+perm(1:end-1))/2 ./ dxF;
    
    dphiCenterGauss = -(permDphiDx(2:end) + permDphiDx(1:end-1))./dxC(2:end-1);

    dphiLeftGauss = permDphiDx(2:end-1) ./ dxC(3:end-1);

    dphiRightGauss = permDphiDx(2:end-1) ./dxC(2:end-2);

    dcCenterGauss = reshape(nA*litersPerCubicMeter*e*vale(1,:) .* repmat(poro(2:end-1), [1, nSpecies]), [nSpecies*(length(dxC)-2), 1]);

    vals(nSpecies*(valsArrayBlockSizePerSpecies)+1:nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies) = ...
        [dphiCenterGauss; ...
        dphiLeftGauss; ...
        dphiRightGauss; ...
        dcCenterGauss];
    
    % Ad hoc modification of vals and indices, in order to prototype higher
    % order scheme. Will make this more official later.
    vals(nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + 1 : nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + nSpecies*2*(N-6)) ...
        = [dPhiLeftLeftAll; dPhiRightRightAll];
    rowInd(nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + 1 : nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + nSpecies*2*(N-6)) ...
        = [rowIndLeftAll, rowIndRightAll];
    colInd(nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + 1 : nSpecies*(valsArrayBlockSizePerSpecies) + 3*N-8 + (N-2)*nSpecies + nSpecies*2*(N-6)) ...
        = [colIndLeftAll, colIndRightAll];

    % Sum of valence times concentration. Length = N.
    gaussEquationSpeciesSum = transpose(sum(reshapedZ(1:nSpecies, :).*vale(1,:)', 1));

    b((nSpecies+1):(nSpecies+1):end) = -nA*litersPerCubicMeter * e * poro(2:end-1) .* gaussEquationSpeciesSum(2:end-1) ...
         - (permDphiDx(2:end) .* (z(3*(nSpecies+1):(nSpecies+1):end) - z(2*(nSpecies+1):(nSpecies+1):end-(nSpecies+1))) ...
         - permDphiDx(1:end-1) .* (z(2*(nSpecies+1):(nSpecies+1):end-(nSpecies+1)) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1))))...
         ./dxC(2:end-1) ...
         - nA*litersPerCubicMeter * e .* backCharge(2:end-1);
 
    %% Solve for the delta and update
    warning('off', 'MATLAB:nearlySingularMatrix')
    if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
        if highPrecision || 1
            A = sparse(rowInd(1:nnz+2*nSpecies*(N-6)), colInd(1:nnz+2*nSpecies*(N-6)), vals(1:nnz+2*nSpecies*(N-6)), nEq, nEq, nnz+2*nSpecies*(N-6));
        else
            sparseMatrixReassign(A, vals, COO2CSCIndices);
        end
        dz = A\b;
    else
        z = NaN(size(z));
        return
    end
    
    % Update the starred variable
    z((nSpecies+1+1):(end-(nSpecies+1))) ...
        = z((nSpecies+1+1):(end-(nSpecies+1))) + dz;
    
    % Exit immediately if first iteration returns imaginary values
    if ~isreal(dz) || any(isnan(dz))
        return
    end
end
%------------- END OF CODE --------------
