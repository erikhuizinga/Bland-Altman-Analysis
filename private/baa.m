function varargout = baa(x, xName, y, yName, a, ...
                         doPlotMD, axMD, ...
                         doPlotMR, axMR, ...
                         MSDType, doPlotMSD1,axMSD1, doPlotMSD2,axMSD2, ...
                         doPlotC, axC, ...
                         doPlotBasicStats, doPlotExtendedStats, ...
                         doPlotRegStats, doConstantRegression, ...
                         doRepeated, assumeCTV)


%% Prepare
if doPlotMD || doPlotC || doPlotMR || doPlotMSD1 || doPlotMSD2
    ax = [axMD; axMR; axMSD1; axMSD2; axC];
    f = get(ax, 'Parent');
    if iscell(f), f = vertcat(f{:}); end
    f = unique(f);  % f can be the handle to one or more figures
else
    f = [];
end

% Check and prepare for repeated measurements
[x, y, n] = prepRep(x, y, doRepeated);


%% Calculate
% Calcultate significance statistics
p = 1 - a/2;
z = Ninv(p);  % Calculate z-statistic for p = 1-alpha/2
t = Tinv(p, n - 1);  % Calculate t-statistics for p with n-1 d.o.f.


% Calculate difference statistics
[muXY, d, varXWithin, varYWithin, loaDCI, loaD, muD, muDCI, eLoaD, ...
    eMuD, sD, polyMuXYD, msePolyMuXYD, sResPolyMuXYD, polyLLoaD, ...
    polyULoaD, m, X, Y] ...
    = statMuS(x, y, 'difference', n, z, t, doConstantRegression, ...
              assumeCTV);


% Calculate mean-difference correlation statistics
[rSMuD, pRSMuD] = corr(muXY, d, 'type', 'Spearman');


% Calculate ratio statistics
if doPlotMR  % This is only calculated if mean-ratio graph is requested
    [~, R, ~, ~, loaRCI, loaR, muR, muRCI, eLoaR, eMuR, sR, ...
        polyMuXYR, msePolyMuXYR, sResPolyMuXYR, polyLLoaR, polyULoaR] ...
        = statMuS(x, y, 'ratio', n, z, t, doConstantRegression, assumeCTV);
    
    % Calculate mean-ratio correlation statistics
    [rSMuR, pRSMuR] = corr(muXY, R, 'type', 'Spearman');
end


% Calculate standard deviation statistics 1
if doPlotMSD1  % This is only calculated is mean-standard deviation graph
               % is requested
    xMSD1 = x;
    
    if strcmp(MSDType, 'separate')
        yMSD1 = x;
    else
        yMSD1 = y;
    end
    
    XYName = [xName, ' and ', yName];
    switch MSDType
        case 'difference'
            xNameMSD1 = ['difference between ', XYName];
        case 'ratio'
            xNameMSD1 = ['ratio between ', XYName];
        otherwise
            xNameMSD1 = xName;
    end
    
    [muMSD1, sMSD1, polyMSD1, msePolyMSD1, polyLLoaMSD1, polyULoaMSD1] ...
        = statMuS(xMSD1, yMSD1, 'SD', z, doConstantRegression, MSDType);
    
    % Calculate mean-standard deviation correlation statistics
    [rMSD1, pRMSD1] = corr(muMSD1, sMSD1, 'type', 'Spearman');
end


% Calculate standard deviation statistics 2
if doPlotMSD2  % This is only calculated when mean-standard deviation graph is requested
    xMSD2 = y;
    yMSD2 = y;
    xNameMSD2 = yName;
    
    [muMSD2, sMSD2, polyMSD2, msePolyMSD2, polyLLoaMSD2, polyULoaMSD2] ...
        = statMuS(xMSD2, yMSD2, 'SD', z, doConstantRegression, MSDType);
    
    % Calculate mean-standard deviation correlation statistics
    [rMSD2, pRMSD2] = corr(muMSD2, sMSD2, 'type', 'Spearman'); %TODO make independent of stats toolbox?
end


% Calculate correlation statistics and linear regression
%TODO linreg for muXY and d
if doPlotC
    [pRhoXY, rhoXY, polyXY, msePXY] = statC(X, Y, z, doConstantRegression);
end


%% Create graphics
% Create correlation plot
if doPlotC
    plotC(axC, X, Y, doPlotBasicStats, pRhoXY, rhoXY, sum(m), xName, yName)
end


% Create mean-difference plot
if doPlotMD
    plotM(axMD, muXY, d, 'difference', 'd', 0, doPlotBasicStats, ...
          loaDCI, pRSMuD, rSMuD, loaD, a, z, muD, muDCI, ...
          doPlotExtendedStats, eLoaD, eMuD, '-', n, xName, yName, ...
          doPlotRegStats, polyMuXYD, msePolyMuXYD, polyLLoaD, ...
          polyULoaD, doConstantRegression, m)
end


% Create mean-ratio plot
if doPlotMR
    plotM(axMR, muXY, R, 'ratio', 'R', 1, doPlotBasicStats, loaRCI, ...
          pRSMuR, rSMuR, loaR, a, z, muR, muRCI, doPlotExtendedStats, ...
          eLoaR,eMuR, '/', n, xName, yName, doPlotRegStats, polyMuXYR, ...
          msePolyMuXYR, polyLLoaR, polyULoaR, doConstantRegression, m)
end


% Create mean-standard deviation plot 1
if doPlotMSD1
    plotM(axMSD1, muMSD1, sMSD1, 'standard deviation', 's', NaN, ...
          doPlotBasicStats, [], pRMSD1,rMSD1, [], [], [], [], [], ...
          doPlotExtendedStats, [], [], 'std', n, xNameMSD1, [], ...
          doPlotRegStats, polyMSD1, msePolyMSD1, polyLLoaMSD1, ...
          polyULoaMSD1, doConstantRegression, m)
end


% Create mean-standard deviation plot 2
if doPlotMSD2
    plotM(axMSD2, muMSD2, sMSD2, 'standard deviation', 's', NaN, ...
          doPlotBasicStats, [], pRMSD2, rMSD2, [], [], [], [], [], ...
          doPlotExtendedStats, [], [], 'std', n, xNameMSD2, [], ...
          doPlotRegStats, polyMSD2, msePolyMSD2, polyLLoaMSD2, ...
          polyULoaMSD2, doConstantRegression, m)
end


%% Set data cursor update function for figure(s)
for f = f(:).'
    dc = datacursormode(f);
    dc.UpdateFcn = @dcUpdateFcn;
    dc.SnapToDataVertex = 'off';
    dc.Enable = 'on';
end


%% Set outputs
% Set difference outputs
out.difference.mu = muD;
out.difference.muCI = muDCI;
out.difference.loa = loaD;
out.difference.loaCI = loaDCI;
out.difference.s = sD;
out.difference.rSMu = rSMuD;
out.difference.pRSMu = pRSMuD;
out.difference.polyMu = polyMuXYD;
out.difference.msePolyMu = msePolyMuXYD;
out.difference.sPolyResidual = sResPolyMuXYD;


% Set ratio outputs
if doPlotMR
    out.ratio.mu = muR;
    out.ratio.muCI = muRCI;
    out.ratio.loa = loaR;
    out.ratio.loaCI = loaRCI;
    out.ratio.s = sR;
    out.ratio.rSMu = rSMuR;
    out.ratio.pRSMu = pRSMuR;
    out.ratio.polyMu = polyMuXYR;
    out.ratio.msePolyMu = msePolyMuXYR;
    out.ratio.sPolyResidual = sResPolyMuXYR;
end


% Set correlation outputs
if doPlotC
    out.correlation.rho = rhoXY;
    out.correlation.p = pRhoXY;
    out.correlation.poly = polyXY;
    out.correlation.polyMSE = msePXY;
end


% Set general outputs
out.x.varWithin = varXWithin;
out.y.varWithin = varYWithin;


% Set final output
varargout = {{out}};
end