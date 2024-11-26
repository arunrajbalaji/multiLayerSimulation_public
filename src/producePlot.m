function producePlot(z, constants, xCenter, xFace, dxC, dxF, speciesToPlot, uniqueSpecies, isElectrodeBC, plotType, figNum)

nSpecies = length(uniqueSpecies);
cMap = linspecer(nSpecies);

if strcmp(plotType, 'concentration')
    figure(figNum)
    clf
    hold on
    for speciesIndex = speciesToPlot
        plot(xCenter, z(speciesIndex:(nSpecies+1):end).*constants.poro(:), '-o', ...
            'Color', cMap(speciesIndex, :), 'MarkerFaceColor', cMap(speciesIndex, :))
    end
    hold off

    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$x$$ (m)', 'interpreter','latex')
    ylabel('$$c$$ (mol/L)', 'interpreter', 'latex')
    legend(uniqueSpecies{speciesToPlot}, 'interpreter', 'latex')

elseif strcmp(plotType, 'waterEq')
    rxnRateCoefficients = constants.rxnRateCoefficients;
    rxnWienCoeff = constants.rxnWienCoeff;
    wienBeta = constants.wienBeta;
    
    figure(figNum)
    clf
    semilogy(xCenter, z(1:(nSpecies+1):end).*z(6:(nSpecies+1):end), '-o', ...
        'Color', cMap(1, :), 'MarkerFaceColor', cMap(1, :))
    dPhidX = (z(3*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1)))./(dxF(1:end-1) + dxF(2:end));
    if size(rxnRateCoefficients, 2) == 8
        rate1 = rxnRateCoefficients(2:end-1,6) ./ rxnRateCoefficients(2:end-1,5);
        rate2 = rxnRateCoefficients(2:end-1,8) ./ rxnRateCoefficients(2:end-1,7);
        finalRate = rate1;
        finalRate(isnan(finalRate)) = rate2(isnan(finalRate));
    elseif size(rxnRateCoefficients, 2) == 6
        finalRate = rxnRateCoefficients(2:end-1,6) ./ rxnRateCoefficients(2:end-1,5);
    end
    
    gammaE = 1./sqrt(exp(rxnWienCoeff(6) * wienBeta * abs(dPhidX)));
    hold on
    semilogy(xCenter(2:end-1), finalRate .* exp(rxnWienCoeff(6) * wienBeta * abs(dPhidX)), ...
        '-o', 'Color', cMap(2, :))
    hold off

    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$x$$ (m)', 'interpreter','latex')
    ylabel('$$K_{\textrm{W}}$$', 'interpreter', 'latex')
    legend({'measured', 'expected'}, 'interpreter', 'latex')
    
    if size(rxnRateCoefficients, 2) == 8
        dissocRate = (rxnRateCoefficients(2:end-1,6) + rxnRateCoefficients(2:end-1,8)) .* exp(rxnWienCoeff(6) * wienBeta * abs(dPhidX));
        recombRate = rxnRateCoefficients(2:end-1,5) + rxnRateCoefficients(2:end-1,7);
    elseif size(rxnRateCoefficients, 2) == 6
        dissocRate = rxnRateCoefficients(2:end-1, 6) .* exp(rxnWienCoeff(6) * wienBeta * abs(dPhidX));
        recombRate = rxnRateCoefficients(2:end-1, 5);
    end
    netRateTotal = (recombRate .* z(nSpecies+1+1:(nSpecies+1):end-(nSpecies+1)).*z(nSpecies+1+6:(nSpecies+1):end-(nSpecies+1)) - dissocRate) .* dxC(2:end-1) .* constants.poro(2:end-1);
    reactionCurrent = sum(netRateTotal(xCenter(2:end-1) > 80e-6 & xCenter(2:end-1) < 120e-6)) * constants.litersPerCubicMeter * constants.nA * constants.e;
    disp(['Reaction current: ', num2str(reactionCurrent), ' A/m2'])

elseif strcmp(plotType, 'charge')
    figure(figNum)
    chargeProfile = constants.backCharge * constants.nA * constants.e * constants.litersPerCubicMeter;
    for speciesIndex = 1:nSpecies
        chargeProfile = chargeProfile + constants.vale(1,speciesIndex) ...
            * z(speciesIndex:(nSpecies+1):end)  * constants.nA * constants.e * constants.litersPerCubicMeter;
    end
    clf
    plot(xCenter, chargeProfile, '-o', 'Color', 'k', 'MarkerFaceColor', 'k')
    hold off

    totalCharge = sum(dxC .* chargeProfile);

    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$x$$ (m)', 'interpreter','latex')
    ylabel('$$\rho$$ (C/m$$^3$$)', 'interpreter', 'latex')
    annotation('textbox', [0.15, 0.10, 0.1, 0.1], 'String', ...
        ['Total Charge: ' num2str(totalCharge) ' C/m$$^2$$'], 'interpreter', 'latex')

elseif strcmp(plotType, 'potential')
    figure(figNum)
    clf
    plot(xCenter(2:end-1), z(2*(nSpecies+1):(nSpecies+1):end-(nSpecies+1)), '-o', ...
        'Color', 'k', 'MarkerFaceColor', 'k')
    hold off

    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$x$$ (m)', 'interpreter','latex')
    ylabel('$$\phi$$ (V)', 'interpreter', 'latex')
    
elseif strcmp(plotType, 'Efield')
    
    dPhidX = -(z(3*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1)))./(dxF(1:end-1) + dxF(2:end));
    
    figure(figNum)
    clf
    plot(xCenter(2:end-1), -dPhidX, '-o', ...
        'Color', 'k', 'MarkerFaceColor', 'k')
    hold off

    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$x$$ (m)', 'interpreter','latex')
    ylabel('$$E$$ (V/m)', 'interpreter', 'latex')
    
elseif strcmp(plotType, 'flux') || strcmp(plotType, 'current')
    % Unpack constants, which are useful for plotting.
    poro = constants.poro;
    tort = constants.tort;
    diff = constants.diff;
    acti = constants.acti;
    vale = constants.vale;
    e = constants.e;
    kb = constants.kb;
    T = constants.T;
    leftElectrodeBC = isElectrodeBC(1);
    rightElectrodeBC = isElectrodeBC(2);
    stericACubed = constants.stericACubed;
    stericOnOffVec = constants.stericOnOffVec;

    totalSpeciesFlux = zeros(length(dxF), nSpecies);
    totalCurrent = zeros(length(dxF), 1);
    
    deltaStern = constants.deltaStern;
    leftExchCurrent = constants.leftIonExchangeCurrent;
    rightExchCurrent = constants.rightIonExchangeCurrent;
    
    faradaicCoeff = deltaStern*e/kb/T;

    wienActivityExponent = constants.wienActivityExponent;

    % Compute summation term for steric effect flux.
    reshapedZ = reshape(z, nSpecies+1, length(dxC));
    stericTermSpeciesSum = transpose(sum(reshapedZ(1:nSpecies, :).*stericOnOffVec, 1));
    rxnWienCoeff = constants.rxnWienCoeff;
    wienBeta = constants.wienBeta;
    for speciesIndex = 1:nSpecies
        % Diffusion flux, length = (N-1), computed at faces (left and
        % right) of all cells that are actively computed.
        diffusionFlux = -1./dxF.*2./(acti(1:end-1, speciesIndex) + acti(2:end, speciesIndex)) ...
            .* (z((nSpecies+1)+speciesIndex:(nSpecies+1):end).*acti(2:end,speciesIndex) - z(speciesIndex:(nSpecies+1):end-(nSpecies+1)).*acti(1:end-1, speciesIndex));
        if leftElectrodeBC
            diffusionFlux(1) = 0;
        end
        if rightElectrodeBC
            diffusionFlux(end) = 0;
        end

        % Electromigration flux, length = (N-1), computed at faces (left
        % and right) of all cells that are actively computed.
        electromigrationFlux = -1/2*(vale(1,speciesIndex)*e/kb/T)*(z(2*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-(nSpecies+1)))./dxF ...
            .* (z(speciesIndex:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+speciesIndex:(nSpecies+1):end));
        if leftElectrodeBC
            electromigrationFlux(1) = 0;
        end
        if rightElectrodeBC
            electromigrationFlux(end) = 0;
        end

        % Steric flux, length = (N-1), computed at faces (left and right)
        % of all cells that are actively computed
        stericFlux = stericOnOffVec(speciesIndex) * 1/2 * log((1 - stericACubed * stericTermSpeciesSum(2:end))./(1 - stericACubed * stericTermSpeciesSum(1:end-1)))./dxF ...
            .* (z(speciesIndex:(nSpecies+1):end-(nSpecies+1)) + z((nSpecies+1)+speciesIndex:(nSpecies+1):end));
        if leftElectrodeBC
            stericFlux(1) = 0;
        end
        if rightElectrodeBC
            stericFlux(end) = 0;
        end

        wienAlphaBeta = rxnWienCoeff(6) * wienBeta;
        wienFieldStrength = wienAlphaBeta .* abs(z(3*(nSpecies+1):(nSpecies+1):end) - z((nSpecies+1):(nSpecies+1):end-2*(nSpecies+1)))./(dxF(2:end) + dxF(1:end-1));
        matrixWienTermDeltaC = wienActivityExponent(speciesIndex) * (1./dxF(2:end-1)) .* (wienFieldStrength(2:end) - wienFieldStrength(1:end-1));
        cInterp = 1/2*(z(speciesIndex:(nSpecies+1):end-(nSpecies+1)) + z(nSpecies+1+speciesIndex:(nSpecies+1):end));
        wienEffectFlux = -matrixWienTermDeltaC .* cInterp(2:end-1);

        totalSpeciesFlux(:, speciesIndex) = (diffusionFlux + electromigrationFlux + stericFlux + [0; wienEffectFlux; 0]) ...
            .* (poro(1:end-1)+poro(2:end))/2.*(diff(1:end-1,speciesIndex)+diff(2:end,speciesIndex))/2./((tort(1:end-1)+tort(2:end))/2)*constants.litersPerCubicMeter;

        % If using electrode boundary condition, use Faradaic reaction term
        % to determine the applied flux at the boundary.
        if leftElectrodeBC
            if vale(1, speciesIndex) ~= 0
                exchangeCurrentConversion = 1/(vale(1,speciesIndex) * constants.e) * 1/constants.nA * 1/constants.litersPerCubicMeter;
                leftFaradaicFlux = -exchangeCurrentConversion*leftExchCurrent(speciesIndex)*2*sinh(faradaicCoeff*(z(2*(nSpecies+1)) - z(nSpecies+1))/dxF(1));
                totalSpeciesFlux(1, speciesIndex) = leftFaradaicFlux * constants.litersPerCubicMeter;
            else
                totalSpeciesFlux(1, speciesIndex) = 0;
            end
        end
        
        if rightElectrodeBC
            if vale(1, speciesIndex) ~= 0
                exchangeCurrentConversion = 1/(vale(1,speciesIndex) * constants.e) * 1/constants.nA * 1/constants.litersPerCubicMeter;
                rightFaradaicFlux = -exchangeCurrentConversion*rightExchCurrent(speciesIndex)*2*sinh(faradaicCoeff*(z(end) - z(end-(nSpecies+1)))/dxF(end));
                totalSpeciesFlux(end, speciesIndex) = rightFaradaicFlux * constants.litersPerCubicMeter;
            else
                totalSpeciesFlux(end, speciesIndex) = 0;
            end
        end
        
        totalCurrent = totalCurrent + vale(1, speciesIndex)*totalSpeciesFlux(:, speciesIndex)*constants.nA*constants.e;
        
    end
    if strcmp(plotType, 'flux')
        figure(figNum)
        clf
        hold on
        for speciesIndex = speciesToPlot
            plot(xFace(3:end-2), totalSpeciesFlux(2:end-1, speciesIndex), '-o', ...
                'Color', cMap(speciesIndex, :), 'MarkerFaceColor', cMap(speciesIndex, :))
        end
        hold off

        set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
        xlabel('$$x$$ (m)', 'interpreter','latex')
        ylabel('$$f$$ (mol/m$$^2\cdot$$s)', 'interpreter', 'latex')
        legend(uniqueSpecies{speciesToPlot}, 'interpreter', 'latex')

    elseif strcmp(plotType, 'current')
        figure(figNum)
        clf
        hold on
        plot(xFace(3:end-2), totalCurrent(2:end-1), '-o', ...
                'Color', 'k', 'MarkerFaceColor', 'k')
        for speciesIndex = speciesToPlot
            plot(xFace(3:end-2), vale(1, speciesIndex)*totalSpeciesFlux(2:end-1, speciesIndex)*constants.nA*constants.e, '-o', ...
                'Color', cMap(speciesIndex, :), 'MarkerFaceColor', cMap(speciesIndex, :))
        end
        hold off

        set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
        xlabel('$$x$$ (m)', 'interpreter','latex')
        ylabel('$$j_x$$ (A/m$$^2$$)', 'interpreter', 'latex')
        legend({'Overall', uniqueSpecies{speciesToPlot}}, 'interpreter', 'latex')
        if median(totalCurrent) ~= 0 && ~isnan(median(totalCurrent))
            ylim([-2*abs(median(totalCurrent)) 2*abs(median(totalCurrent))])
        end
        
%         disp(['Boundary current computed for (K+, HSO4-, SO42-): (' ...
%             num2str(vale(1, 4)*totalSpeciesFlux(2, 4)*constants.nA*constants.e) ', ' ...
%             num2str(vale(1, 3)*totalSpeciesFlux(end-1, 3)*constants.nA*constants.e) ', ' ...
%             num2str(vale(1, 6)*totalSpeciesFlux(end-1, 6)*constants.nA*constants.e) ')'])
        
        dPhidX = (z(2*(nSpecies+1):(nSpecies+1):end)-z((nSpecies+1):(nSpecies+1):end-(nSpecies+1)))./dxF;
        leftEDLIndex = find(dPhidX<max(dPhidX)/1e3, 1, 'first');
        rightEDLIndex = length(dPhidX)-find(flipud(dPhidX)<max(dPhidX)/1e3, 1, 'first');
%         leftEDLIndex = find(xCenter>=2e-9, 1, 'first');
%         rightEDLIndex = length(xCenter)-find(flipud(xCenter)<=8e-2+180e-6-2e-9, 1, 'first');
        
        disp(['Computed potential drop, without electrode EDL: ', ...
            num2str(z(rightEDLIndex*(nSpecies+1)) - z(leftEDLIndex*(nSpecies+1))), ' V'])
        
        disp(['Median current: ', num2str(median(totalCurrent)), ' A/m2.'])
        disp(['H+ current: ', num2str(median(vale(1, 1)*totalSpeciesFlux(1:ceil(length(xCenter)/2), 1)*constants.nA*constants.e)), ' A/m2.'])
        disp(['OH-: ', num2str(median(vale(1, 6)*totalSpeciesFlux(ceil(length(xCenter)/2):end, 6)*constants.nA*constants.e)), ' A/m2.'])
        disp(['Maximum E-field: ', num2str(max(abs(dPhidX))), ' V/m.'])
    end
else
    disp('Error: specified plot type is not supported')
end
end