function plotM(axMY, ...
    x,y, ...
    sMYLongName,sMYName, ...
    y0, ...
    doPlotBasicStats, ...
    loaCI,pRSMuY,rSMuY,loa,a,z,muY,muYCI, ...
    doPlotExtendedStats,eLoa,eMuY, ...
    strYFun,n,xName,yName, ...
    doPlotRegStats,polyXY,msePXY,polyLLoa,polyULoa,doConReg, ...
    m ...
    )
% mean-statistic plot, Y referring to the statistic


doMSD = strcmp(sMYLongName,'standard deviation');
if doMSD && ~isscalar(axMY)
    % plot 1
    plotM(axMY(1), ...
        x(:,1),y(:,1), ...
        sMYLongName,sMYName, ...
        y0, ...
        doPlotBasicStats, ...
        loaCI,pRSMuY(1),rSMuY(1),loa,a,z,muY,muYCI, ...
        doPlotExtendedStats,eLoa,eMuY, ...
        strYFun,n,xName,xName, ... % note the double xName
        doPlotRegStats, ...
        polyXY(1,:),msePXY(1),polyLLoa(1,:),polyULoa(1,:), ...
        doConReg, ...
        m)
    % plot 2
    plotM(axMY(2), ...
        x(:,2),y(:,2), ...
        sMYLongName,sMYName, ...
        y0, ...
        doPlotBasicStats, ...
        loaCI,pRSMuY(2),rSMuY(2),loa,a,z,muY,muYCI, ...
        doPlotExtendedStats,eLoa,eMuY, ...
        strYFun,n,yName,yName, ... % idem for yName
        doPlotRegStats, ...
        polyXY(2,:),msePXY(2),polyLLoa(2,:),polyULoa(2,:), ...
        doConReg, ...
        m)
    return
end

% preparation
axes(axMY)
hold(axMY,'on')
legEntries = gobjects(0);
strP = num2str(100*(1-a));
M = sum(m); % number of observation pairs

% mean-Y plot
sM = scatter(x,y);
sM.ZData = 1:n;
UserData = dcStruct([],'µ',sMYName,'i',[],@dcXYZ); %TODO check µ
if any(m>1)
    % repeated measurements
    sM.UserData{1} = UserData;
    sM.UserData{2} = m;
else
    % no repeated measurements
    sM.UserData = UserData;
end
xl = xlim;
yl = ylim;
padding = range(yl)/20; % distance to enhance visibility

% line of equality
if ~doMSD
    line0 = refline(0,y0);
    line0.Color = [.75 .75 .75];
    line0.LineStyle = '--';
    line0.DisplayName = 'line of equality';
    line0.UserData = dcStruct('M1 = M2',[],[],[],[],@dcXYZ);
    legEntries(end+1) = line0;
    
    % adjust y limits
    yl = [ ...
        min(yl(1), y0-padding), ...
        max(yl(2), y0+padding) ...
        ];
    ylim(yl)
end

% plot statistics
if doPlotBasicStats
    % add correlation to legend
    if pRSMuY<1e-4
        % as recommended by BMJ 1996;312:572
        % http://dx.doi.org/10.1136/bmj.312.7030.572
        strPRSMuD = sprintf('\\itp\\rm < 0.0001');
    else
        strPRSMuD = sprintf('\\itp\\rm = %.2f',pRSMuY);
    end
    sM.DisplayName = sprintf( ...
        '\\itr_s\\rm = %.2f (%s)',rSMuY,strPRSMuD);
    legEntries(end+1) = sM;
    
    if ~doPlotRegStats && ~doMSD
        % adjust y limits
        yl = ylim;
        padding = range(yl)/20; % distance to enhance visibility
        yl = [ ...
            min(yl(1), loa(1)-padding), ...
            max(yl(2), loa(2)+padding) ...
            ];
        ylim(yl)
        
        % lower LOA line
        lineLLoa = refline(0,loa(1));
        lineLLoa.Color = 'k';
        lineLLoa.UserData = dcStruct([],[], ...
            sprintf('lower %s%% LOA',strP), ...
            [], ...
            [strP '% CI = [' num2str(loaCI(1,1)) ', ' ...
            num2str(loaCI(2,1)) ']'], ...
            @dcXYZ);
        text(xl(2)-padding/2,loa(1)+padding/2, ...
            sprintf( ...
            '$\\overline{%s}-%.2fs_%s$',sMYName,z,sMYName ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % upper limit of agreement line
        lineULoa = refline(0,loa(2));
        lineULoa.Color = 'k';
        lineULoa.UserData = dcStruct([],[], ...
            sprintf('upper %s%% LOA',strP), ...
            [], ...
            [strP '% CI = [' ...
            num2str(loaCI(1,2)) ', ' num2str(loaCI(2,2)) ']'], ...
            @dcXYZ);
        text(xl(2)-padding/2,loa(2)+padding/2, ...
            sprintf( ...
            '$\\overline{%s}+%.2fs_%s$',sMYName,z,sMYName ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % mean Y line
        muYLine = refline(0,muY);
        muYLine.Color = 'k';
        text(xl(2)-padding/2,muY+padding/2, ...
            sprintf('$\\overline{%s}$',sMYName), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        muYLine.UserData = dcStruct([],[],['mean ' sMYLongName], [], ...
            [strP '% CI = [' ...
            num2str(muYCI(1)) ', ' num2str(muYCI(2)) ']'], ...
            @dcXYZ);
    end
    
    % plot additional statistics
    if doPlotExtendedStats && ~doMSD
        % lower LOA errorbar
        lE = errorbar(xl(2),loa(1),eLoa,'r');
        lE.UserData = lineLLoa.UserData;
        
        % mean Y errorbar
        muYE = errorbar(xl(2),muY,eMuY,'r');
        muYE.UserData = muYLine.UserData;
        
        % upper LOA errorbar
        uE = errorbar(xl(2),loa(2),eLoa,'r');
        uE.UserData = lineULoa.UserData;
        
        % adjust y limits
        yl = ylim;
        padding = range(yl)/20; % distance to enhance visibility
        yl = [ ...
            min(yl(1), loaCI(1,1)-padding), ...
            max(yl(2), loaCI(2,2)+padding) ...
            ];
        ylim(yl)
    end
end

% plot regression statistics
if doPlotRegStats
    % polyXY: polynomial coefficients
    % msePXY: MSE of polynomial fit
    % polyLLoa: polynomial coefficients of the lower LOA
    % polyULoa: idem for the upper LOA
    
    % mean statistic regression line
    yRegXY = polyval(polyXY,xlim);
    regLine = line(xlim,yRegXY);
    regLine.Color = 'k';
    regLine.DisplayName = ...
        sprintf('regression with %s%% LOA', ...
        strP);
    regLine.UserData = dcStruct( ...
        sprintf('f_bias(µ) = %s×µ + %s', ... %TODO check × and µ
        num2str(polyXY(1)),num2str(polyXY(2))), ...
        'µ',sMYName,[],[], ... %TODO check µ
        @dcXYZ);
    legEntries(end+1) = regLine;
    
    % mean statistic regression text
    text(xl(2)-padding/2,yRegXY(2)+sign(polyXY(1))*padding/2, ...
        sprintf('$f_%s(\\mu)$',sMYName), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','right')
    
    % lower LOA line
    yLLoa = polyval(polyLLoa,xlim);
    regLineLLoa = line(xlim,yLLoa);
    regLineLLoa.Color = 'k';
    if doConReg
        strLoaFun = sprintf('s_{res}');
        strLoaDC = sprintf('%.2f×s_res',z); %TODO check ×
    else
        strLoaFun = sprintf('\\sqrt{\\pi/2}f_{res}(\\mu)');
        strLoaDC = sprintf('%.2f×f_res(µ)',z*sqrt(pi/2)); %TODO check symbols
    end
    regLineLLoa.UserData = dcStruct( ...
        sprintf( ...
        'f_lower(µ) - %.2f×f_res(µ) = %s×µ + %s', ... %TODO check symbols
        sqrt(pi/2)*z, ...
        num2str(polyLLoa(1)), num2str(polyLLoa(2)) ...
        ), ...
        'µ', sprintf('lower %s%% LOA',strP), [], [], ... %TODO check µ
        @dcXYZ);
    
    % upper LOA line
    yULoa = polyval(polyULoa,xlim);
    regLineULoa = line(xlim,yULoa);
    regLineULoa.Color = 'k';
    regLineULoa.UserData = dcStruct( ...
        sprintf( ...
        'f_upper(µ) + %s = %s×µ + %s', ... %TODO check symbols
        strLoaDC, ...
        num2str(polyULoa(1)),num2str(polyULoa(2)) ...
        ), ...
        'µ', sprintf('upper %s%% LOA',strP), [], [], ... %TODO check µ
        @dcXYZ);
    
    % adjust y limits
    yl = ylim;
    padding = range(yl)/20; % distance to enhance visibility
    yl = [ ...
        min(yl(1), min(yLLoa)-padding), ...
        max(yl(2), max(yULoa)+padding) ...
        ];
    ylim(yl)
    
    % lower LOA text
    text(xl(2)-padding/2,yLLoa(2)+sign(polyLLoa(1))*padding/2, ...
        sprintf( ...
        '$f_%s(\\mu)-%.2f%s$',sMYName,z,strLoaFun ...
        ), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','right')
    
    % upper LOA text
    text(xl(2)-padding/2,yULoa(2)+sign(polyULoa(1))*padding/2, ...
        sprintf( ...
        '$f_%s(\\mu)+%.2f%s$',sMYName,z,strLoaFun ...
        ), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','right')
end

% axes labels
if n==M
    strSubObs = sprintf('\\itn = \\Sigmam\\bf = %u subjects and observation pairs',n);
elseif doMSD
    strSubObs = sprintf('\\itn\\bf = %u subjects, %u observations',n,M);
else
    strSubObs = sprintf('\\itn\\bf = %u subjects, %u observation pairs',n,M);
end
if doMSD
    xlabel('subject mean')
    ylabel(sMYLongName)
    title( ...
        sprintf( ...
        'Mean-%s plot of (%s):\n\\rm%s', ...
        sMYLongName,strSubObs,xName ...
        ) ...
        )
else
    xlabel(sprintf('mean \\itµ = (M_1+M_2)/2')) %TODO check µ
    strYLabel = sprintf('%s \\it%s = M_1%sM_2',sMYLongName,sMYName,strYFun);
    ylabel(strYLabel)
    title(sprintf( ...
        ['Mean-%s plot of (%s):\n' ...
        ' \\rm\\itM_1\\rm: %s, \\itM_2\\rm: %s'], ...
        sMYLongName,strSubObs,xName,yName))
end

% legend
if ~isempty(legEntries)
    legend(legEntries,'Location','SouthWest')
end

% reorder plot children
axMY.Children = axMY.Children([end 1:end-1]);

% % equalise axes
% axis equal

% set outer position of axes
setOP(axMY)
end
