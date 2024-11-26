function [constants] = genPhysConstArrays(uniqueSpecies, layerInfo, constants, rxnInfo, xCenter, interfaceLocs)
% genPhysConstArrays: Function generates arrays that correspond to the 1D
% domain grid, for physical constants that vary in space due to variation
% between layers. Arrays are 2D if they are also species-dependent.
% Formatting constants in this way helps with code vectorization.
%
% [constants] = genPhysConstArrays(uniqueSpecies, layerInfo, constants, xCenter)
%
% Inputs:
%       uniqueSpecies   - ordered list of unique species names, from
%                         parsing function.
%       layerInfo       - layer info struct, as output by parsing function.
%       constants       - constants struct, as output by parsing function.
%                         Arrays will be stored here and then returned as a
%                         larger struct.
%       xCenter         - List of cell centers, where values are stored.
%                         Should not include an interface between layers.
%       interfaceLocs   - List of interfaces between contiguous layers.
%                         Must be of size one more than number of layers.
%
% Outputs:
%       constants.poro            - 1D array of porosity
%       constants.tort            - 1D array of tortuosity
%       constants.perm            - 1D array of electric permittivity
%       constants.backCharge      - 1D array of background charge density
%       constants.diff            - 2D array of species diffusivities
%       constants.acti            - 2D array of species tortuosities
%       constants.vale            - 2D array of species diffusivities
%       constants.initVal         - 2D array of species initial values
%       constants.stericACubed    - Shortcut parameter for steric effect
%       constants.wienBeta        - Appears in second wien Effect
%       constants.leftNoFlux      - Logical array, is species no-flux at
%                                   left
%       constants.rightNoFlux     - Logical array, is species no-flux at
%                                   right
%       constants.stericOnOffVec  - Logical array, is steric effect active
%                                   for species or not. Set based on
%                                   non-zero valence.
%       constants.rxnRateCoefficients   - Array of reaction rate coefficients
%       constants.rxnInputOutputIndices - Indices of reactants and products
%       constants.rxnWienCoeff          - wienCoefficient (appears inside
%                                         exponential, multiples entire
%                                         argument).
%
%       NOTE: only new outputs in addition to the existing memebrs of the
%       constants struct are listed above. For all 1D arrays, primary
%       dimension corresponds to grid. For all 2d arrays, rows correspond
%       to grid, and collumns correspond to different species. See
%       parseInputFile.m
%
% Other m-files required: none
% MAT-files required: none
%
% See also: parseInputFile.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 29-August-2021
%------------- BEGIN CODE --------------
% Unpack preliminary data
nSpecies = length(uniqueSpecies);
nLayers = length(interfaceLocs) - 1;

% Non-species-dependent properties (porosity, tortuosity, permeability
constants.poro = zeros(size(xCenter));
constants.tort = zeros(size(xCenter));
constants.perm = zeros(size(xCenter));
constants.backCharge = zeros(size(xCenter));

for ii = 1:nLayers
%    constants.poro((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
%        = layerInfo(ii).poro;
%    constants.tort((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
%        = layerInfo(ii).tort;
%    constants.perm((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
%        = layerInfo(ii).perm;
%     constants.backCharge((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
%         = layerInfo(ii).backCharge;
    
   if nLayers == 1
       constants.poro(:) = layerInfo(ii).poro;
       constants.tort(:) = layerInfo(ii).tort;
       constants.perm(:) = layerInfo(ii).perm;
       constants.backCharge(:) = layerInfo(ii).backCharge;
   elseif ii == 1
       constants.poro = constants.poro ...
           + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).poro;
       constants.tort = constants.tort ...
           + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).tort;
       constants.perm = constants.perm ...
           + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).perm;
       constants.backCharge = constants.backCharge ...
           + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).backCharge;
   elseif ii == nLayers
       constants.poro = constants.poro ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).poro;
       constants.tort = constants.tort ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).tort;
       constants.perm = constants.perm ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).perm;
       constants.backCharge = constants.backCharge ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
           *layerInfo(ii).backCharge;
   else
       constants.poro = constants.poro ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
           - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
           *layerInfo(ii).poro;
       constants.tort = constants.tort ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
           - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
           *layerInfo(ii).tort;
       constants.perm = constants.perm ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
           - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
           *layerInfo(ii).perm;
       constants.backCharge = constants.backCharge ...
           + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
           - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
           *layerInfo(ii).backCharge;
   end
end

% Species-dependent properties (diffusivity and activity, valence and
% initialization value)
constants.diff = zeros(length(xCenter), nSpecies);
constants.acti = zeros(length(xCenter), nSpecies);
constants.vale = zeros(length(xCenter), nSpecies);
constants.initVal = zeros(length(xCenter), nSpecies);

for ii = 1:nLayers
    for jj = 1:nSpecies
        % If species is not present, then diffusivity should be zero
        xPointsInLayerIndices = (xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1));
        constants.diff(xPointsInLayerIndices,jj)...
            = layerInfo(ii).diff(jj);
        
        % We apply a hyperbolic tangent smoothing function to the activity
        % coefficient as a function of space, so that gradients of the
        % activity coefficient are well defined (not mesh dependent)
        % If species is not present, then activity should be one (in LOG)
        if nLayers == 1
            constants.acti(:,jj) = layerInfo(ii).acti(jj);
        elseif ii == 1
            constants.acti(:,jj) = constants.acti(:,jj) ...
                + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
                *layerInfo(ii).acti(jj);
        elseif ii == nLayers
            constants.acti(:,jj) = constants.acti(:,jj) ...
                + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
                *layerInfo(ii).acti(jj);
        else
            constants.acti(:,jj) = constants.acti(:,jj) ...
                + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
                - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
                *layerInfo(ii).acti(jj);
        end
        
        % If species is not present, then valence should be zero (but
        % doesn't matter, since diffusivity is already set to zero). Here,
        % will just set all values for all of space equal to the value for
        % the last time valence is specified. Valence should be specified
        % as same for each species in every layer.
        if layerInfo(ii).vale(jj) ~= 0
            constants.vale(:,jj) = layerInfo(ii).vale(jj);
        end
        
        % constants.initVal((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1)),jj)...
        %    = layerInfo(ii).init(jj);
         if nLayers == 1
              constants.initVal(:, jj) = layerInfo(ii).init(jj);
         elseif ii == 1
              constants.initVal(:, jj) = constants.initVal(:, jj) ...
                 + (-1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) + 1/2)...
                 * layerInfo(ii).init(jj);
         elseif ii == nLayers
             constants.initVal(:, jj) = constants.initVal(:, jj) ...
                 + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2)...
                 *layerInfo(ii).init(jj);
         else
             constants.initVal(:, jj) = constants.initVal(:, jj) ...
                 + (1/2*tanh((xCenter - interfaceLocs(ii))/constants.deltaActivity) + 1/2 ...
                 - 1/2*tanh((xCenter - interfaceLocs(ii+1))/constants.deltaActivity) - 1/2)...
                 *layerInfo(ii).init(jj);
         end
    end
end

% Steric Effects, Second Wien Effect parameters.
constants.stericACubed = constants.litersPerCubicMeter * constants.nA * constants.stericA^3;
constants.stericOnOffVec = [constants.vale(1,:).*sign(constants.vale(1,:))./max(ones(size(constants.vale(1,:))), abs(constants.vale(1,:)))]';

l_B = constants.e^2/(4*pi*constants.waterPerm*constants.kb*constants.T);
constants.wienBeta = l_B*constants.e/(constants.kb*constants.T);
% ASSUMPTION: NO VARIATION IN SPACE OR BETWEEN LAYERS, WE SIMPLY TAKE THE
% VALUES SPECIFIED IN THE FIRST LAYER. VALUES IN OTHER LAYERS ARE NEGLECTED.
constants.wienActivityExponent = layerInfo(1).wienActivityExponent;

% Reaction handling
% First, we determine the overall set of unique reactions for entire domain.
% A reaction is considered a match for another reaction if they have the
% same reactants, products, reaction rate coefficient, and Wien
% coefficient.
uniqueReactionSet = [];
for ii = 1:nLayers
    for jj = 1:size(rxnInfo(ii).speciesRxns, 1)
        uniqueReactionSet = [uniqueReactionSet; rxnInfo(ii).speciesRxns(jj, :)];
    end
end
uniqueReactionSet = unique(uniqueReactionSet, 'stable', 'rows');
nReactionsTotal = size(uniqueReactionSet, 1);

% Reaction rate coefficients are stored in vectorized form, to allow for
% vectorized computations. It is possible for these values to vary in
% space. Input/output indices and wien coefficients are stored separately.
constants.rxnRateCoefficients = zeros(length(xCenter), nReactionsTotal);
constants.rxnInputOutputIndices = zeros(nReactionsTotal, 4);
constants.rxnWienCoeff = zeros(nReactionsTotal, 1);
% Loop over all grid points
for ii = 1:length(xCenter)
    % Determien which layer we are in
    layerNumber = max(find((xCenter(ii) > interfaceLocs)));
    for jj = 1:size(rxnInfo(layerNumber).speciesRxns, 1)
        % Determine a perfect match. In case of same reaction specified
        % multiple times in a single layer, all the reaction information
        % will be collapsed into a single reaction.
        matchingReactionIndexInUniqueReactionSet = find(sum(uniqueReactionSet == rxnInfo(layerNumber).speciesRxns(jj,:), 2) == 6);
        
        % Construct array of reaction rate coefficients
        constants.rxnRateCoefficients(ii, matchingReactionIndexInUniqueReactionSet) = rxnInfo(layerNumber).speciesRxns(jj, 5);
        
        % These parameters do not vary in x, but we compute them and set
        % them here (redundantly, since this in a loop over grid points).
        % Cost is minimal, since this code is only every run once.
        constants.rxnInputOutputIndices(matchingReactionIndexInUniqueReactionSet, :) = rxnInfo(layerNumber).speciesRxns(jj, 1:4);
        constants.rxnWienCoeff(matchingReactionIndexInUniqueReactionSet) = rxnInfo(layerNumber).speciesRxns(jj, 6);
    end
end

%------------- END OF CODE --------------
