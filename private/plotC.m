function scatterC = plotC(axC, x, y, ...
    doPlotBasicStats, pRhoXY, rhoXY, n, ...
    xName, yName, ...
    scatterName, m, titleObsStr)
% Create correlation plot of observations not accounting for repeated
% measurements


% Calculate total number of observations
N = sum(m);


% Prepare graphics objects
axes(axC)
legEntries = gobjects(0);


% Determine if honeycomb plot is requested
isHoneycomb = strcmpi(scatterName, 'honeycomb');


% Scatter y against x
scatterFunction = getScatterFunction(scatterName);
scatterC = scatterFunction(x, y);
if ~isHoneycomb
    % Set ZData to the subject number
    scatterC.ZData = getZData(x, n, N, m);
    
    % Set data for custom data tip
    scatterC.UserData = dcStruct([], 'M1', 'M2', 'j', [], @dcXYZ);
end


% Add statistics
if doPlotBasicStats
    % Add correlation to legend
    if pRhoXY < 1e-4
        % Change string when pRhoXY < 1e-4, as recommended by BMJ
        % 1996;312:572: http://dx.doi.org/10.1136/bmj.312.7030.572
        strPRhoXY = sprintf('\\itp\\rm < 0.0001');
        
    else
        strPRhoXY = sprintf('\\itp\\rm = %.2f', pRhoXY);
    end
    
    scatterC.DisplayName = sprintf('\\rho = %.2f (%s)', rhoXY, strPRhoXY);
    legEntries(end + 1) = scatterC;
end


% Add y = x reference line
eqLine = refline(1, 0);  % Line at 45°, because axis are equal
eqLine.Color = [getLineColor(scatterName), .5];
eqLine.LineStyle = '--';
eqLine.DisplayName = 'line of equality';
eqLine.UserData = dcStruct([], [], [], [], 'M2 = M1', @dcXYZ);
legEntries(end + 1) = eqLine;


% Add axes labels
xlabel('M_1')
ylabel('M_2')
scatterName = lower(scatterName);
scatterName(1) = upper(scatterName(1));
if n == N
    titleSuffix = sprintf( ...
        '\\itn\\rm\\bf = %u subjects and observation pairs', n);
    
else
    titleSuffix = sprintf(['\\itn\\rm\\bf = %u subjects, ' ...
                           '\\Sigma\\itm\\rm\\bf = %u observation ' ...
                           'pairs'], n, N);
end
title(sprintf([scatterName, ' plot ', titleObsStr, ' (%s)\n' ...
               '\\rm\\itM_1\\rm: %s, \\itM_2\\rm: %s'], ...
               titleSuffix, xName, yName))

% Add legend
legend(legEntries, 'Location', 'SouthEast')

% Reorder plot children
if ~isHoneycomb
    axC.Children = axC.Children([end, 1 : end-1]);
end

% Equalise axes
axis tight
axis equal

% Set outer position
setOP(axC)
end
