function [leftBCValues, rightBCValues] = computeExternalSolution(z, dxC, dxF, constants, originalLeftBC, originalRightBC, lResLeft, lResRight)

% Useful constants
rightSaltC = sum(abs(constants.vale(end,:)) .* originalRightBC');
rightAmbiD = sum(constants.diff(end,:) .* abs(constants.vale(end,:)) .* originalRightBC') / rightSaltC;
leftSaltC = sum(abs(constants.vale(1,:)) .* originalLeftBC');
leftAmbiD = sum(constants.diff(1,:) .* abs(constants.vale(1,:)) .* originalLeftBC') / leftSaltC;
nSpecies = length(originalLeftBC);

% Constructing total salt concentration
totalSalt = zeros(length(dxC), 1);
for sIndex = 1:nSpecies
    totalSalt = totalSalt + abs(constants.vale(1,sIndex)) * z(sIndex:(nSpecies+1):end);
end

% Measuring fluxes at the boundaries
zRightFlux = -rightAmbiD * median((totalSalt(end-9:end)-totalSalt(end-10:end-1))./dxF(end-9:end));
zLeftFlux = -leftAmbiD * median((totalSalt(2:11)-totalSalt(1:10))./dxF(1:10));

calcRightC = rightSaltC + min(lResRight * zRightFlux/rightAmbiD, 0);
rightBCFactor = calcRightC/rightSaltC;
rightBCValues = rightBCFactor * originalRightBC;

calcLeftC = leftSaltC - max(lResLeft * zLeftFlux/leftAmbiD, 0);
leftBCFactor =  calcLeftC/leftSaltC;
leftBCValues = leftBCFactor * originalLeftBC;

end