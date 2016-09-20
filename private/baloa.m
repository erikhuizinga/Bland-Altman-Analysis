function varargout = baloa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotMR, axMR, ...
    doPlotMSD, axMSD, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtendedStats, ...
    doPlotRegStats, doConstantRegression, ...
    doRepeated)
%% preparation
if doPlotMD || doPlotC || doPlotMR || doPlotMSD
    % if doMDPlot, f(1) = axMD.Parent; end
    % if doCPlot, f(2) = axC.Parent; end
    ax = [axMD;axMR;axMSD;axC];
    f = get(ax,'Parent');
    if iscell(f), f = vertcat(f{:}); end
    f = unique(f);
    % f can be the handle to one or more figures
else
    f = [];
end

% check and prepare for repeated measurements
[x,y,n] = prepRep(x,y,doRepeated);

%% calculations
% significance statistics
p = 1-a/2;
z = Ninv(p); % inverse normal distribution at p = 1-alpha/2
t = Tinv(p,n-1); % inverse t-distribution at p

% difference statistics
[muXY,d,varXWithin,varYWithin,loaDCI,loaD,muD,muDCI,eLoaD,eMuD,sD, ...
    polyMuXYD,msePolyMuXYD,sResPolyMuXYD,polyLLoaD,polyULoaD,m,X,Y] = ...
    statMuS(x,y,'difference',n,z,t,doConstantRegression);

% mean-difference correlation statistics
[rSMuD,pRSMuD] = corr(muXY,d,'type','Spearman'); %TODO make independent of stats toolbox?

% ratio statistics
if doPlotMR % only calculated when mean-ratio graph is requested
    [~,R,~,~,loaRCI,loaR,muR,muRCI,eLoaR,eMuR,sR, ...
        polyMuXYR,msePolyMuXYR,sResPolyMuXYR,polyLLoaR,polyULoaR] = ...
        statMuS(x,y,'ratio',n,z,t,doConstantRegression);
    
    % mean-ratio correlation statistics
    [rSMuR,pRSMuR] = corr(muXY,R,'type','Spearman'); %TODO make independent of stats toolbox?
end

% standard deviation statistics
if doPlotMSD % only calculated when mean-standard deviation graph is requested
    [muX,muY,sX,sY, ...
        polyMSDX,msePolyMSDX,~,polyLLoaMSDX,polyULoaMSDX, ...
        polyMSDY,msePolyMSDY,~,polyLLoaMSDY,polyULoaMSDY ...
        ] = statMuS(x,y,'SD',z,doConstantRegression);
    % concatenate
    muMSD = [muX,muY];
    sMSD = [sX,sY];
    polyMSD = [polyMSDX;polyMSDY];
    msePolyMSD = [msePolyMSDX,msePolyMSDY];
    % sResPolyMSD = [sResPolyMSDX,sResPolyMSDY]; %TODO not used?
    polyLLoaMSD = [polyLLoaMSDX;polyLLoaMSDY];
    polyULoaMSD = [polyULoaMSDX;polyULoaMSDY];
    
    % mean-standard deviation correlation statistics
    [rSMuSX,pRSMuSX] = corr(muX,sX,'type','Spearman'); %TODO make independent of stats toolbox?
    [rSMuSY,pRSMuSY] = corr(muY,sY,'type','Spearman'); %TODO make independent of stats toolbox?
    rMSD = [rSMuSX,rSMuSY];
    pRMSD = [pRSMuSX,pRSMuSY];
end

% correlation statistics and linear regression %TODO linreg for muXY and d
if doPlotC
    [pRhoXY,rhoXY,polyXY,msePXY] = statC(X,Y,z,doConstantRegression);
end

%% graphics
% correlation plot
if doPlotC
    plotC(axC,X,Y,doPlotBasicStats,pRhoXY,rhoXY,sum(m),xName,yName)
end

% mean-difference plot
if doPlotMD
    % plotMD(axMD,muXY,d,doRatio,doPlotBasicStats,loaCI,pRSMuD,rSMuD,loa,a,z,muD,muDCI,doPlotExtStats,eLoa,eMuD,doPlotLS,n,xName,yName)
    plotM(axMD,muXY,d,'difference','d',0,doPlotBasicStats,loaDCI, ...
        pRSMuD,rSMuD,loaD,a,z,muD,muDCI,doPlotExtendedStats,eLoaD,eMuD, ...
        '-',n,xName,yName, ...
        doPlotRegStats,polyMuXYD,msePolyMuXYD,polyLLoaD,polyULoaD,doConstantRegression, ...
        m)
end

% mean-ratio plot
if doPlotMR
    plotM(axMR,muXY,R,'ratio','R',1,doPlotBasicStats,loaRCI,pRSMuR, ...
        rSMuR,loaR,a,z,muR,muRCI,doPlotExtendedStats,eLoaR,eMuR, ...
        '/',n,xName,yName, ...
        doPlotRegStats,polyMuXYR,msePolyMuXYR,polyLLoaR,polyULoaR,doConstantRegression, ...
        m)
end

% mean-standard deviation plot
if doPlotMSD
    plotM(axMSD, muMSD,sMSD, 'standard deviation','s', ...
        NaN, doPlotBasicStats, ...
        [], ...
        pRMSD,rMSD, ...
        [],[],[],[],[],doPlotExtendedStats,[],[], ...
        'std',n,xName,yName, ...
        doPlotRegStats, ...
        polyMSD,msePolyMSD,polyLLoaMSD,polyULoaMSD,doConstantRegression, ...
        m)
end

%% set data cursor update function for figure(s)
for f = f(:).'
    dc = datacursormode(f);
    dc.UpdateFcn = @dcUpdateFcn;
    dc.SnapToDataVertex = 'off';
    dc.Enable = 'on';
end

%% output
% difference outputs
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

% ratio outputs
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

% correlation outputs
if doPlotC
    out.correlation.rho = rhoXY;
    out.correlation.p = pRhoXY;
    out.correlation.poly = polyXY;
    out.correlation.polyMSE = msePXY;
end

% general outputs
out.x.varWithin = varXWithin;
out.y.varWithin = varYWithin;

% final output
varargout = {{out}};
end
