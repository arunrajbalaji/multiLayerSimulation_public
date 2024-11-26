function [rowInd, colInd, nnz, nEq, valsArrayBlockSizePerSpecies, COO2CSCIndices, A] = makeSparseMatrixIndexVectors_fourthOrder(nSpecies, N, nFourthOrderCells)
% makeSparseMatrixIndexVectors: Produces the row and column indec vectors
% for the matrix that must be inverted in the linear solve, for each
% iteration within each time step.
%
% Examine rowInd and coldInd vectors to see the order in which matrix
% coefficient values should be added to the 'vals' vector in doTimeStep.m. 
% Make sure that order of entries matches, so that array is constructed
% properly.
%
%   [rowInd, colInd, nnz, nEq] = makeSparseMatrixIndexVectors(nSpecies, N, constants)
%
% Inputs:
%
%       nSpecies        - Number of uique species
%       N               - Number of cell centers
%
% Outputs:
%       rowInd          - Matrix row indices, in correct order
%       colInd          - Matrix column indices, in correct order
%       nnz             - Number of non-zeros in the matrix
%       nEq             - Number of  unknowns (number of equations)
%       valsArrayBlockSizePerSpecies
%                       - Number of sparse entries per species. Useful for
%                       construting the corresponding 'vals' matrix in
%                       doTimeStep.m
%
%
% Other m-files required: none
% MAT-files required: none
%
% See also: doTimeStep.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2022
%------------- BEGIN CODE --------------

% Data allocation

% For a given species, we expect dependence on the concentration of each
% individual species and the electric potential at the center, left, and
% right cells.
valsArrayBlockSizePerSpecies = (nSpecies+1)*(3*N-8);

% Number of cells for which a higher order scheme in space must be used,
% requiring values from one more cell over on either side. All
% concentrations and potential in the two cells over are required.
valsArrayFourthOrderPerSpecies = (nSpecies+1)*(2 * nFourthOrderCells);
first4thOrderInd = ceil(N/2) - ceil(nFourthOrderCells/2);
last4thOrderInd = first4thOrderInd + nFourthOrderCells - 1;

% The total number of non-zeros is the number of species multiplied by the
% amount per species (comptued above), plus the tri-diagonal system for
% Poisson's equation, plus the supplementary fourth-order cells
nnz = nSpecies*valsArrayBlockSizePerSpecies + (3*N-8) + nSpecies*(N-2) + nSpecies * valsArrayFourthOrderPerSpecies;

totalNumberPotential = (3*N-8) + nSpecies*(N-2);

nEq = (N-2)*(nSpecies + 1);     % Edge cells not computed, Dirichlet BC
rowInd = zeros(nnz,1);
colInd = zeros(nnz,1);

% Loop over iterations
%% Determine LHS matrix A
% Loop over species
for jj = 1:nSpecies
    rowInd((jj-1)*(valsArrayBlockSizePerSpecies)+1:jj*(valsArrayBlockSizePerSpecies)) = ...
        [repmat((jj:(nSpecies+1):(N-2)*(nSpecies+1))', [nSpecies, 1]); ...
        repmat(((nSpecies+1)+jj:(nSpecies+1):(N-2)*(nSpecies+1))', [nSpecies, 1]); ...
        repmat((jj:(nSpecies+1):(N-3)*(nSpecies+1))', [nSpecies, 1]); ...
        (jj:(nSpecies+1):(N-2)*(nSpecies+1))'; ...
        ((nSpecies+1)+jj:(nSpecies+1):(N-2)*(nSpecies+1))'; ...
        (jj:(nSpecies+1):(N-3)*(nSpecies+1))'];
    
    colInd((jj-1)*(valsArrayBlockSizePerSpecies)+1:jj*(valsArrayBlockSizePerSpecies)) = ...
        [reshape((repmat((1:nSpecies)', [1, N-2]) + repmat(1:(nSpecies+1):((nSpecies+1)*(N-2)), [nSpecies, 1]) - 1)', [nSpecies*(N-2), 1]); ...
        reshape((repmat((1:nSpecies)', [1, N-3]) + repmat(1:(nSpecies+1):((nSpecies+1)*(N-3)), [nSpecies, 1]) - 1)', [nSpecies*(N-3), 1]); ...
        reshape((repmat((1:nSpecies)', [1, N-3]) + repmat((nSpecies+1)+1:(nSpecies+1):((nSpecies+1)*(N-2)), [nSpecies, 1]) - 1)', [nSpecies*(N-3), 1]); ...
        ((nSpecies+1):(nSpecies+1):(N-2)*(nSpecies+1))'; ...
        ((nSpecies+1):(nSpecies+1):(N-3)*(nSpecies+1))'; ...
        (2*(nSpecies+1):(nSpecies+1):(N-2)*(nSpecies+1))'];
end

% Add rows/columns for Poisson's equation for electric potential
rowInd(nSpecies*(valsArrayBlockSizePerSpecies)+1:nSpecies*(valsArrayBlockSizePerSpecies)+totalNumberPotential) = ...
    [(nSpecies+1:nSpecies+1:(N-2)*(nSpecies+1))'; ...
    (2*(nSpecies+1):nSpecies+1:(N-2)*(nSpecies+1))'; ...
    (nSpecies+1:nSpecies+1:(N-3)*(nSpecies+1))'; ...
    repmat((nSpecies+1:nSpecies+1:(N-2)*(nSpecies+1))', [nSpecies, 1])];
colInd(nSpecies*(valsArrayBlockSizePerSpecies)+1:nSpecies*(valsArrayBlockSizePerSpecies)+totalNumberPotential) = ...
    [(nSpecies+1:nSpecies+1:(N-2)*(nSpecies+1))'; ...
    ((nSpecies+1):nSpecies+1:(N-3)*(nSpecies+1))'; ...
    (2*(nSpecies+1):nSpecies+1:(N-2)*(nSpecies+1))'; ...
    reshape((repmat((1:nSpecies)', [1, N-2]) + repmat(1:(nSpecies+1):((nSpecies+1)*(N-2)), [nSpecies, 1]) - 1)', [nSpecies*(N-2), 1])];

entriesSoFar = nSpecies*(valsArrayBlockSizePerSpecies)+totalNumberPotential;

% Add rows/columns for fourth order treatment of specific cells
for jj = 1:nSpecies
    rowInd(entriesSoFar + (jj-1)*valsArrayFourthOrderPerSpecies + 1:entriesSoFar + (jj)*valsArrayFourthOrderPerSpecies) = ...
        [repmat(((nSpecies+1)*(first4thOrderInd-1) + jj : (nSpecies+1) : last4thOrderInd*(nSpecies+1))', [nSpecies,1]); ...
        repmat(((nSpecies+1)*(first4thOrderInd-1) + jj : (nSpecies+1) : last4thOrderInd*(nSpecies+1))', [nSpecies,1]); ...
        ((nSpecies+1)*(first4thOrderInd-1) + jj : (nSpecies+1) : last4thOrderInd*(nSpecies+1))'; ...
        ((nSpecies+1)*(first4thOrderInd-1) + jj : (nSpecies+1) : last4thOrderInd*(nSpecies+1))'];
    
    colInd(entriesSoFar + (jj-1)*valsArrayFourthOrderPerSpecies + 1:entriesSoFar + (jj)*valsArrayFourthOrderPerSpecies) = ...
        [reshape((repmat((1:nSpecies)', [1, nFourthOrderCells]) + repmat((first4thOrderInd-1 - 2)*(nSpecies+1)+1:(nSpecies+1):((nSpecies+1)*(last4thOrderInd - 2)), [nSpecies, 1]) - 1)', [nSpecies*(nFourthOrderCells), 1]); ...
        reshape((repmat((1:nSpecies)', [1, nFourthOrderCells]) + repmat((first4thOrderInd-1 + 2)*(nSpecies+1)+1:(nSpecies+1):((nSpecies+1)*(last4thOrderInd + 2)), [nSpecies, 1]) - 1)', [nSpecies*(nFourthOrderCells), 1]); ...
        ((nSpecies+1)*(first4thOrderInd) : (nSpecies+1) : last4thOrderInd*(nSpecies+1))'; ...
        ((nSpecies+1)*(first4thOrderInd) : (nSpecies+1) : last4thOrderInd*(nSpecies+1))'];
end


% Initialize space in memory.
A = sparse(rowInd, colInd, ones(size(rowInd)), nEq, nEq, nnz);

% Compute conversion vector for going from COO format to CSC format. This is
% used in the sparseMatrixReassign mex file, which can quickly update
% sparse matrix memory values in place. Performed in a slow and inefficient
% way, but only needs to be done once. Could be computed analytically, if
% desired...
COO2CSCIndices = int32(zeros(nnz, 1));
for nnzIndex = 1:nnz

    numCellsToLeft = (colInd(nnzIndex)-1 - mod(colInd(nnzIndex)-1, nSpecies+1))/(nSpecies+1);
    numEntriesInLowerColumns = max(numCellsToLeft-1, 0)*((nSpecies+1)^2 + 2*(nSpecies*(nSpecies+1) + 1)) ...
        + (numCellsToLeft > 0)*((nSpecies+1)^2 + nSpecies*(nSpecies+1) + 1);

    numLowerColumnsInSameCell = colInd(nnzIndex)-1 - numCellsToLeft*(nSpecies+1);
    numCellsAbove = (rowInd(nnzIndex)-1 - mod(rowInd(nnzIndex)-1, nSpecies+1))/(nSpecies+1);
    numCellsBelow = N-2 - numCellsAbove - 1;
    numCellsToRight = N-2 - numCellsToLeft - 1;
    numEntriesInLowerColumnsInSameCell = numLowerColumnsInSameCell * ...
        (nSpecies+1 + nSpecies*(numCellsToLeft > 0) + nSpecies*(numCellsToRight > 0));
    
    variableIdentifier = numLowerColumnsInSameCell;
    equationIdentifier = rowInd(nnzIndex)-1 - numCellsAbove*(nSpecies+1);

    if numCellsAbove < numCellsToLeft % check if located in upper band
            numEntriesInSameColumn = equationIdentifier;
    elseif numCellsAbove == numCellsToLeft % check if located in middle band
        if numCellsAbove > 0 % Make sure not entry (1,1)
            if variableIdentifier == nSpecies % check if potential
                numEntriesInSameColumn = equationIdentifier + nSpecies+1;
            else
                numEntriesInSameColumn = equationIdentifier + nSpecies;
            end
        else
            numEntriesInSameColumn = equationIdentifier;
        end
    elseif numCellsAbove > numCellsToLeft % check if located in lower band
        if numCellsAbove > 1 % Make sure not entry (1,2)
            if variableIdentifier == nSpecies % check if potential
                numEntriesInSameColumn = equationIdentifier + 2*(nSpecies+1);
            else
                numEntriesInSameColumn = equationIdentifier + 2*nSpecies + 1;
            end
        else
            if variableIdentifier == nSpecies % check if potential
                numEntriesInSameColumn = equationIdentifier + (nSpecies+1);
            else
                numEntriesInSameColumn = equationIdentifier + nSpecies + 1;
            end
        end
    end

    COO2CSCIndices(nnzIndex) = numEntriesInLowerColumns + numEntriesInLowerColumnsInSameCell ...
        + numEntriesInSameColumn + 1;
end

end

%------------- END OF CODE --------------