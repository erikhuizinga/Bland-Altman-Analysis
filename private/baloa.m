function varargout = baloa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotMR, axMR, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtStats, ...
    doPlotRegStats, doConReg, ...
    doPlotLS, ...
    doRepeated)
%% preparation
if doPlotMD || doPlotC || doPlotMR
    % if doMDPlot, f(1) = axMD.Parent; end
    % if doCPlot, f(2) = axC.Parent; end
    ax = [axMD;axMR;axC];
    f = get(ax,'Parent');
    if iscell(f), f = vertcat(f{:}); end
    f = unique(f);
    % f can be the handle to one or more figures
    % if numel(f)>1, multiFig = true; end
else
    f = [];
end

if doRepeated % BAA for repeated measurements
    % distinguish the number of replicates
    if iscell(x) || iscell(y)
        % Unequal number of replicates. This is certain , because if
        % original input was a cell and contained equal number of
        % observations per subject, it was converted to a matrix.
        %TODO
    else % x and y are not cells, thus matrices
        if size(x,2) == size(y,2)
            % equal number of replicates
            
        else
            % unequal number of replicates
            %TODO
        end
    end
end

% keep only values that can be used in calculations
lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
n = nnz(lok);
xok = x(lok);
yok = y(lok);

%% calculations
% mean statistics
muXY = mean([xok,yok],2);

% significance statistics
p = 1-a/2;
z = Ninv(p); % inverse normal distribution at p = 1-alpha/2
t = Tinv(p,n-1); % inverse t-distribution at p

% difference statistics
d = xok-yok;
[loaDCI,loaD,muD,muDCI,eLoaD,eMuD,sD,polyMuXYD,msePolyMuXYD, ...
    sResPolyMuXYD,polyLLoaD,polyULoaD] = statY(muXY,d,n,z,t,doConReg);

% ratio statistics
R = xok./yok;
[loaRCI,loaR,muR,muRCI,eLoaR,eMuR,sR,polyMuXYR,msePolyMuXYR, ...
    sResPolyMuXYR,polyLLoaR,polyULoaR] = statY(muXY,R,n,z,t,doConReg);

% mean-difference correlation statistics
[rSMuD,pRSMuD] = corr(muXY,d,'type','Spearman'); %TODO make independent of stats toolbox?

% mean-ratio correlation statistics
[rSMuR,pRSMuR] = corr(muXY,R,'type','Spearman'); %TODO make independent of stats toolbox?

% correlation statistics and linear regression %TODO linreg for muXY and d
[pRhoXY,rhoXY,polyXY,msePXY] = statC(xok,yok,z,doConReg);

%% graphics
% correlation plot
if doPlotC
    plotC(axC,xok,yok,doPlotBasicStats,pRhoXY,rhoXY,doPlotLS,polyXY, ...
        msePXY,n,xName,yName)
end

% mean-difference plot
if doPlotMD
    % plotMD(axMD,muXY,d,doRatio,doPlotBasicStats,loaCI,pRSMuD,rSMuD,loa,a,z,muD,muDCI,doPlotExtStats,eLoa,eMuD,doPlotLS,n,xName,yName)
    plotM(axMD,muXY,d,'difference','d',0,doPlotBasicStats,loaDCI, ...
        pRSMuD,rSMuD,loaD,a,z,muD,muDCI,doPlotExtStats,eLoaD,eMuD, ...
        doPlotLS,'-',n,xName,yName, ...
        doPlotRegStats,polyMuXYD,msePolyMuXYD,polyLLoaD,polyULoaD,doConReg)
end

% mean-ratio plot
if doPlotMR
    plotM(axMR,muXY,R,'ratio','R',1,doPlotBasicStats,loaRCI,pRSMuR, ...
        rSMuR,loaR,a,z,muR,muRCI,doPlotExtStats,eLoaR,eMuR,doPlotLS, ...
        '/',n,xName,yName, ...
        doPlotRegStats,polyMuXYR,msePolyMuXYR,polyLLoaR,polyULoaR,doConReg)
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

% final output
varargout = {{out}};
end