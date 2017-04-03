function scatterM = plotM( ...
               axMY, x, y, sMYLongName, sMYName, y0, doPlotBasicStats, ...
               loaCI, pRSMuY, rSMuY, loa, a, z, muY, muYCI, ...
               doPlotExtendedStats, eLoa, eMuY, strYFun, n, xName, ...
               yName, doPlotRegStats, polyXY, msePXY, polyLLoa, ...
               polyULoa, doConReg, m, scatterName)
% Create mean-statistic plot, y referring to the statistic

% Determine if mean-SD plot
doMSD = strcmp(sMYLongName, 'standard deviation');

% Prepare some stuff
axes(axMY)
hold(axMY, 'on')
legEntries = gobjects(0);
strP = num2str(100 * (1 - a));
N = sum(m); % number of observation pairs


% Create mean-y plot
scatterFunction = getScatterFunction(scatterName);
scatterM = scatterFunction(x, y);
if strcmpi(scatterName, 'scatter')
    scatterM.ZData = 1 : n;
    UserData = dcStruct([], 'µ', sMYName, 'i', [], @dcXYZ); %TODO check µ
    
else
    bar = colorbar;
    ylabel(bar, 'counts', 'Rotation', -90, 'VerticalAlignment', 'bottom')
    UserData = struct.empty;
end

if any(m > 1)  % Handle repeated measurements
    scatterM.UserData{1} = UserData;
    scatterM.UserData{2} = m;
    
else  % Handle no repeated measurements
    scatterM.UserData = UserData;
end

xl = xlim;
yl = ylim;
padding = range(yl)/20;  % Add a distance to enhance visibility


% Add line of equality
if ~doMSD
    line0 = refline(0, y0);
    line0.Color = [.75, .75, .75];
    line0.LineStyle = '--';
    line0.DisplayName = 'line of equality';
    line0.UserData = dcStruct('M1 = M2', [], [], [], [], @dcXYZ);
    legEntries(end + 1) = line0;
    
    % Adjust y limits
    yl = [min(yl(1), y0 - padding), ...
          max(yl(2), y0 + padding)];
    ylim(yl)
end

% Plot statistics
if doPlotBasicStats
    % Add correlation to legend
    if pRSMuY < 1e-4
        % Change string when pRhoXY < 1e-4, as recommended by BMJ
        % 1996;312:572: http://dx.doi.org/10.1136/bmj.312.7030.572
        strPRSMuD = sprintf('\\itp\\rm < 0.0001');
    else
        strPRSMuD = sprintf('\\itp\\rm = %.2g', pRSMuY);
    end
    scatterM.DisplayName = sprintf('\\itr_s\\rm = %.2f (%s)', rSMuY, strPRSMuD);
    legEntries(end + 1) = scatterM;
    
    if ~doPlotRegStats && ~doMSD
        % Adjust y limits
        yl = ylim;
        padding = range(yl) / 20;  % Add a distance to enhance visibility
        yl = [min(yl(1), loa(1) - padding), ...
              max(yl(2), loa(2) + padding)];
        ylim(yl)
        
        % Add lower LOA line
        lineLLoa = refline(0, loa(1));
        lineLLoa.Color = 'k';
        lineLLoa.UserData = dcStruct([], [], ...
                                     sprintf('lower %s%% LOA', strP), ...
                                     [], [strP, '% CI = [', ...
                                          num2str(loaCI(1, 1)), ', ', ...
                                          num2str(loaCI(2, 1)),  ']'], ...
                                     @dcXYZ);
        text(xl(2) - padding/2, loa(1) + padding/2, ...
             sprintf('$\\overline{%s}-%.2fs_%s=%s$', ...
                     sMYName, z, sMYName, num2str(loa(1))), ...
             'Interpreter', 'latex', 'HorizontalAlignment', 'right')
        
        % Add upper LOA line
        lineULoa = refline(0, loa(2));
        lineULoa.Color = 'k';
        lineULoa.UserData = dcStruct([], [], ...
                                     sprintf('upper %s%% LOA', strP), ...
                                     [], [strP, '% CI = [', ...
                                          num2str(loaCI(1,2)), ', ', ...
                                          num2str(loaCI(2,2)), ']'], ...
                                     @dcXYZ);
        text(xl(2) - padding/2, loa(2) + padding/2, ...
             sprintf('$\\overline{%s}+%.2fs_%s=%s$', ...
                     sMYName, z, sMYName, num2str(loa(2))), ...
             'Interpreter', 'latex', 'HorizontalAlignment', 'right')
        
        % Add mean Y line
        muYLine = refline(0, muY);
        muYLine.Color = 'k';
        text(xl(2) - padding/2, muY + padding/2, ...
             sprintf('$\\overline{%s}=%s$', sMYName, num2str(muY)), ...
             'Interpreter', 'latex', 'HorizontalAlignment','right')
        muYLine.UserData = dcStruct([], [], ['mean ', sMYLongName], [], ...
                                    [strP, '% CI = [', ...
                                     num2str(muYCI(1)), ', ', ...
                                     num2str(muYCI(2)), ']'], @dcXYZ);
    end
    
    
    % Plot additional statistics
    if doPlotExtendedStats && ~doMSD
        % Add lower LOA errorbar
        lE = errorbar(xl(2), loa(1), eLoa, 'r');
        lE.UserData = lineLLoa.UserData;
        
        % Add mean Y errorbar
        muYE = errorbar(xl(2), muY, eMuY, 'r');
        muYE.UserData = muYLine.UserData;
        
        % Add upper LOA errorbar
        uE = errorbar(xl(2), loa(2), eLoa, 'r');
        uE.UserData = lineULoa.UserData;
        
        % Adjust y limits
        yl = ylim;
        padding = range(yl) / 20;  % Add a distance to enhance visibility
        yl = [min(yl(1), loaCI(1, 1) - padding), ...
              max(yl(2), loaCI(2, 2) + padding)];
        ylim(yl)
    end
end


% Plot regression statistics
if doPlotRegStats
    % polyXY: polynomial coefficients
    % msePXY: MSE of polynomial fit
    % polyLLoa: polynomial coefficients of the lower LOA
    % polyULoa: idem for the upper LOA
    
    % Add mean statistic regression line
    yRegXY = polyval(polyXY, xlim);
    regLine = line(xlim, yRegXY);
    regLine.Color = 'k';
    regLine.DisplayName = sprintf('regression with %s%% LOA', strP);
    regLine.UserData = dcStruct( ...
        sprintf('f_bias(µ) = %s×µ + %s', ...  %TODO check × (times) and µ
        num2str(polyXY(1)), num2str(polyXY(2))), 'µ', ...  %TODO check µ
        sMYName, [], [], @dcXYZ);
    legEntries(end + 1) = regLine;
    
    % Add mean statistic regression text
    text(xl(2) - padding/2, yRegXY(2) + sign(polyXY(1)) * padding/2, ...
        sprintf('$f_{\\overline{%s}}(\\mu)$', sMYName), 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right')
    
    % Add lower LOA line
    yLLoa = polyval(polyLLoa, xlim);
    regLineLLoa = line(xlim, yLLoa);
    regLineLLoa.Color = 'k';
    if doConReg
        strLoaFun = sprintf('s_{\\varepsilon}');
        
    else
        strLoaFun = sprintf('\\sqrt{\\pi/2}f_{\\varepsilon}(\\mu)');
    end
    regLineLLoa.UserData = dcStruct( ...
        sprintf('f_lower(µ) = %s×µ + %s', ...  %TODO check × (times) and µ
        num2str(polyLLoa(1)), num2str(polyLLoa(2))), 'µ', ...  %TODO check µ
        sprintf('lower %s%% LOA', strP), [], [], @dcXYZ);
    
    % Add upper LOA line
    yULoa = polyval(polyULoa, xlim);
    regLineULoa = line(xlim, yULoa);
    regLineULoa.Color = 'k';
    regLineULoa.UserData = dcStruct( ...
        sprintf('f_upper(µ) = %s×µ + %s', ...  %TODO check × (times) and µ
        num2str(polyULoa(1)), num2str(polyULoa(2))), ...
        'µ', ...  %TODO check µ
        sprintf('upper %s%% LOA', strP), [], [], @dcXYZ);
    
    % Adjust y limits
    yl = ylim;
    padding = range(yl) / 20;  % Add a distance to enhance visibility
    yl = [min(yl(1), min(yLLoa) - padding), ...
          max(yl(2), max(yULoa) + padding)];
    ylim(yl)
    
    % Add lower LOA text
    text(xl(2) - padding/2, yLLoa(2) + sign(polyLLoa(1)) * padding/2, ...
         sprintf('$f_{\\overline{%s}}(\\mu)-%.2f%s$', sMYName, z, ...
                 strLoaFun), ...
         'Interpreter', 'latex', 'HorizontalAlignment', 'right')
    
    % Add upper LOA text
    text(xl(2) - padding/2, yULoa(2) + sign(polyULoa(1)) * padding/2, ...
         sprintf('$f_{\\overline{%s}}(\\mu)+%.2f%s$', sMYName, z, ...
                 strLoaFun), ...
         'Interpreter', 'latex', 'HorizontalAlignment', 'right')
end

% Add axes labels
if n == N
    strSubObs = sprintf( ...
        '\\itn\\bf = %u subjects and observation pairs', n);
    
elseif doMSD
    strSubObs = sprintf( ...
        '\\itn\\bf = %u subjects, \\it\\Sigmam\\bf = %u observations', ...
        n, N);
    
else
    strSubObs = sprintf( ...
        '\\itn\\bf = %u subjects, \\it\\Sigmam\\bf = %u observation pairs', ...
        n, N);
end

if doMSD
    xlabel('subject mean')
    ylabel(sMYLongName)
    title(sprintf('Mean-%s plot of (%s):\n\\rm%s', sMYLongName, ...
                  strSubObs, xName))
    
else
    xlabel(sprintf('mean \\itµ = (M_1+M_2)/2'))  %TODO check µ
    strYLabel = sprintf('%s \\it%s = M_1%sM_2', sMYLongName, sMYName, ...
                        strYFun);
    ylabel(strYLabel)
    title(sprintf(['Mean-%s plot of (%s):\n' ...
                   ' \\rm\\itM_1\\rm: %s, \\itM_2\\rm: %s'], ...
                   sMYLongName, strSubObs, xName, yName))
end

% Add legend
if ~isempty(legEntries)
    legend(legEntries, 'Location', 'SouthWest')
end

% Reorder plot children
axMY.Children = axMY.Children([end, 1 : end-1]);

% Set outer position of axes
setOP(axMY)
end
