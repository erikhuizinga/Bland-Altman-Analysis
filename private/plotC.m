function plotC(axC, x, y, doPlotBasicStats, pRhoXY, rhoXY, n, xName, yName)
% Create correlation plot of observations not accounting for repeated
% measurements

% Prepare
axes(axC)
legEntries = gobjects(0);

% Plot y against x
sC = scatter(x, y);
sC.ZData = 1 : n;
sC.UserData = dcStruct([], 'M1', 'M2', 'j', [], @dcXYZ);

if doPlotBasicStats
    % Add correlation to legend
    if pRhoXY < 1e-4
        % Change string when pRhoXY < 1e-4, as recommended by BMJ
        % 1996;312:572: http://dx.doi.org/10.1136/bmj.312.7030.572
        strPRhoXY = sprintf('\\itp\\rm < 0.0001');
        
    else
        strPRhoXY = sprintf('\\itp\\rm = %.2f', pRhoXY);
    end
    
    sC.DisplayName = sprintf('\\rho = %.2f (%s)', rhoXY, strPRhoXY);
    legEntries(end + 1) = sC;
end

% Add y = x reference line
eqLine = refline(1, 0);  % Line at 45°, because axis are equal
eqLine.Color = [.75, .75, .75];
eqLine.LineStyle = '--';
eqLine.DisplayName = 'line of equality';
eqLine.UserData = dcStruct([], [], [], [], 'M2 = M1', @dcXYZ);
legEntries(end + 1) = eqLine;

% Add axes labels
xlabel('M_1')
ylabel('M_2')
title(sprintf(['Scatter plot of (%u observation pairs):\n' ...
               ' \\rm\\itM_1\\rm: %s, \\itM_2\\rm: %s'], n, xName, yName))

% Add legend
legend(legEntries, 'Location', 'SouthEast')

% Reorder plot children
axC.Children = axC.Children([end, 1 : end-1]);

% Equalise axes
axis tight
axis equal

% Set outer position
setOP(axC)
end
