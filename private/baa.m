function varargout = baa(x, xName, y, yName, a, ...
                         doPlotMD, axMD, doPlotMR, axMR, MSDType, ...
                         doPlotMSD1, axMSD1, doPlotMSD2, axMSD2, ...
                         doPlotC, axC, doPlotBasicStats, ...
                         doPlotExtendedStats, doPlotRegStats, ...
                         doConstantRegression, doRepeated, assumeCTV, ...
                         scatterName)


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
[muXY, D, varXWithin, varYWithin, loaDCI, loaD, muD, muDCI, eLoaD, ...
    eMuD, sD, polyMuXYD, msePolyMuXYD, sResPolyMuXYD, polyLLoaD, ...
    polyULoaD, m, X, Y] ...
    = statMuS(x, y, 'difference', n, z, t, doConstantRegression, ...
              assumeCTV);


% Calculate mean-difference correlation statistics
[rSMuD, pRSMuD] = corr(muXY, D, 'type', 'Spearman');


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
    scatterC = plotC( ...
        axC, X, Y, ...
        doPlotBasicStats, pRhoXY, rhoXY, sum(m), ...
        xName, yName, ...
        scatterName);
end


% Create mean-difference plot
if doPlotMD
    scatterMD = plotM( ...
          axMD, muXY, D, 'difference', 'D', 0, doPlotBasicStats, ...
          loaDCI, pRSMuD, rSMuD, loaD, a, z, muD, muDCI, ...
          doPlotExtendedStats, eLoaD, eMuD, '-', n, xName, yName, ...
          doPlotRegStats, polyMuXYD, msePolyMuXYD, polyLLoaD, ...
          polyULoaD, doConstantRegression, m, scatterName);
end


% Create mean-ratio plot
if doPlotMR
    scatterMR = plotM( ...
          axMR, muXY, R, 'ratio', 'R', 1, doPlotBasicStats, loaRCI, ...
          pRSMuR, rSMuR, loaR, a, z, muR, muRCI, doPlotExtendedStats, ...
          eLoaR,eMuR, '/', n, xName, yName, doPlotRegStats, polyMuXYR, ...
          msePolyMuXYR, polyLLoaR, polyULoaR, doConstantRegression, m, ...
          scatterName);
end


% Create mean-standard deviation plot 1
if doPlotMSD1
    % Determine string of SD
    switch MSDType
        case 'difference'
            SDString1 = 's_D';
            
        case 'ratio'
            SDString1 = 's_R';
            
        otherwise
            SDString1 = 's';
    end
    
    % Create mean-SD plot 1
    scatterMSD1 = plotM( ...
          axMSD1, muMSD1, sMSD1, 'standard deviation', SDString1, NaN, ...
          doPlotBasicStats, [], pRMSD1,rMSD1, [], [], [], [], [], ...
          doPlotExtendedStats, [], [], 'std', n, xNameMSD1, [], ...
          doPlotRegStats, polyMSD1, msePolyMSD1, polyLLoaMSD1, ...
          polyULoaMSD1, doConstantRegression, m, scatterName);
end


% Create mean-standard deviation plot 2
if doPlotMSD2
    % Determine string of SD
    switch MSDType
        case 'difference'
            SDString2 = 's_D';
            
        case 'ratio'
            SDString2 = 's_R';
            
        otherwise
            SDString2 = 's';
    end
    
    % Create mean-SD plot 2
    scatterMSD2 = plotM( ...
          axMSD2, muMSD2, sMSD2, 'standard deviation', SDString2, NaN, ...
          doPlotBasicStats, [], pRMSD2, rMSD2, [], [], [], [], [], ...
          doPlotExtendedStats, [], [], 'std', n, xNameMSD2, [], ...
          doPlotRegStats, polyMSD2, msePolyMSD2, polyLLoaMSD2, ...
          polyULoaMSD2, doConstantRegression, m, scatterName);
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
out.difference.bias = muD;  % bias (mean difference (D))
out.difference.biasCI = muDCI;  % confidence interval (CI) of bias
out.difference.loa = loaD;  % limits of agreement (LOA) of D
out.difference.loaCI = loaDCI;  % CI of the LOA of D
out.difference.std = sD;  % standard deviation (SD) of D
out.difference.D = D;  % difference to plot against mean
out.difference.Spearman.r = rSMuD;  % Spearman rank correlation of D and mean (mu)
out.difference.Spearman.p = pRSMuD;  % p-value of rSMuD
out.difference.poly.bias = polyMuXYD;  % simple linear regression of D on mu
out.difference.poly.mse = msePolyMuXYD;  % mean squared error (MSE) of polyMuXYD
out.difference.poly.stde = sResPolyMuXYD;  % SD of the residuals
out.difference.poly.loa = [polyLLoaD.', polyULoaD.'];  % simple linear regression of LOA on mu


% Set ratio outputs
if doPlotMR
    out.ratio.bias = muR;  % mean ratio
    out.ratio.biasCI = muRCI;
    out.ratio.loa = loaR;
    out.ratio.loaCI = loaRCI;
    out.ratio.std = sR;
    out.ratio.R = R;  % ratio to plot against mean
    out.ratio.Spearman.r = rSMuR;
    out.ratio.Spearman.p = pRSMuR;
    out.ratio.poly.bias = polyMuXYR;
    out.ratio.poly.mse = msePolyMuXYR;
    out.ratio.poly.stde = sResPolyMuXYR;
    out.ratio.poly.loa = [polyLLoaR.', polyULoaR.'];
    out.graphics.ratio.handle = scatterMR;  % handle to scatter object
end


% Set correlation outputs
if doPlotC
    out.xy.Pearson.rho = rhoXY;  % Pearson correlation of x with y
    out.xy.Pearson.p = pRhoXY;  % p-value of rhoXY
    out.xy.poly.xy = polyXY;  % simple linear regression of y on x
    out.xy.poly.mse = msePXY;  % MSE of the regression
    out.graphics.xy.handle = scatterC;  % handle to scatter object
end


% Set general outputs
out.xy.x = x;  % x-values used in the calculations
out.xy.y = y;  % y-values used in the calculations
out.xy.mu = muXY;  % mean values to plot against
out.n = n;  % number of subjects


% Set repeated measurements outputs
if doRepeated
    out.xy.varw.x = varXWithin;  % within-subject variance for x
    out.xy.varw.y = varYWithin;  % within-subject variance for y
    out.m = m;  % number of observations per subject
    out.N = sum(m);  % total number of observations
end


% Set remaining graphic outputs
if doPlotMD
    out.grahpics.difference.handle = scatterMD;
end
if doPlotMSD1
    if doPlotMSD2
        out.grahpics.std.handle = [scatterMSD1; scatterMSD2];
    else
        out.graphics.std.handle = scatterMSD1;
    end
end

% Set final output
varargout = {{out}};
end