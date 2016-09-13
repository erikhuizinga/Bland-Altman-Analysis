function plotC(axC,xok,yok,doPlotBasicStats,pRhoXY,rhoXY,doPlotLS,polyXY,msePXY,n,xName,yName)
% correlation plot

% preparation
axes(axC)
legEntries = gobjects(0);

% plot y against x
sC = scatter(xok,yok);
sC.UserData = dcStruct([],'M1','M2',[],@dcXY);

if doPlotBasicStats
    % add correlation to legend
    if pRhoXY<1e-4
        % as recommended by BMJ 1996;312:572
        % http://dx.doi.org/10.1136/bmj.312.7030.572
        strPRhoXY = sprintf('\\itp\\rm < 0.0001');
    else
        strPRhoXY = sprintf('\\itp\\rm = %.2f',pRhoXY);
    end
    sC.DisplayName = sprintf('\\rho = %.2f (%s)',rhoXY,strPRhoXY);
    legEntries(end+1) = sC;
end

% plot least-squares line
if doPlotLS
    XYLine = refline; %TODO make independent of stats toolbox?
    XYLine.Color = 'r';
    legEntries(end+1) = XYLine;
    XYLine.DisplayName = 'least-squares';
    XYLine.UserData = dcStruct( ...
        ['M2 = ' num2str(polyXY(1)) 'M1 + ' num2str(polyXY(2))], ...
        [], [], ...
        ['MSE = ' num2str(msePXY)], ...
        @dcXY);
end

% y = x reference line
eqLine = refline(1,0); % 45° degree, because of axis equal
eqLine.Color = [.75 .75 .75];
eqLine.LineStyle = '--';
eqLine.DisplayName = 'line of equality';
eqLine.UserData = dcStruct([],[],[],'M2 = M1',@dcXY);
legEntries(end+1) = eqLine;

% axes labels
xlabel('M_1')
ylabel('M_2')
title(sprintf(['Scatter plot of (%u observations):\n' ...
    ' \\rm\\itM_1\\rm: %s\n \\itM_2\\rm: %s'],n,xName,yName))

% legend
legend(legEntries,'Location','SouthEast')

% reorder plot children
axC.Children = axC.Children([end 1:end-1]);

% equalise axes
axis tight
axis equal
end