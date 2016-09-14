function plotM(axMY, ...
    x,y, ...
    sMYLongName,sMYName, ...
    y0, ...
    doPlotBasicStats, ...
    loaCI,pRSMuY,rSMuY,loa,a,z,muY,muYCI, ...
    doPlotExtStats,eLoa,eMuY, ...
    doPlotLS, ...
    strYFun,n,xName,yName, ...
    doPlotRegStats,polyXY,msePXY,polyLLoa,polyULoa)
% mean-statistic plot, Y referring to the statistic

% preparation
axes(axMY)
hold(axMY,'on')
legEntries = gobjects(0);
strP = num2str(100*(1-a));

% mean-Y plot
sM = scatter(x,y);
sM.UserData = dcStruct([],'�',sMYName,[],@dcXY); %TODO check �
xl = xlim;

% line of equality
line0 = refline(0,y0);
line0.Color = [.75 .75 .75];
line0.LineStyle = '--';
line0.DisplayName = 'line of equality';
line0.UserData = dcStruct('M1 = M2',[],[],[],@dcXY);
legEntries(end+1) = line0;

% adjust y limits
yl = ylim;
padding = range(yl)/20; % distance to enhance visibility
yl = [ ...
    min(yl(1), y0-padding), ...
    max(yl(2), y0+padding) ...
    ];
ylim(yl)

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
    
    if ~doPlotRegStats
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
            [strP '% CI = [' num2str(loaCI(1,1)) ', ' ...
            num2str(loaCI(2,1)) ']'], ...
            @dcXY);
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
            [strP '% CI = [' ...
            num2str(loaCI(1,2)) ', ' num2str(loaCI(2,2)) ']'], ...
            @dcXY);
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
        muYLine.UserData = dcStruct([],[],['mean ' sMYLongName], ...
            [strP '% CI = [' ...
            num2str(muYCI(1)) ', ' num2str(muYCI(2)) ']'], ...
            @dcXY);
    end
    
    % plot additional statistics
    if doPlotExtStats
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
        sprintf('f_%s(�) = %s׵ + %s', ... %TODO check � and �
        sMYName,num2str(polyXY(1)),num2str(polyXY(2))), ...
        '�',sMYName, ... %TODO check �
        'simple linear regression of difference on mean', ...
        @dcXY);
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
    regLineLLoa.UserData = dcStruct( ...
        sprintf('lower %s%% LOA',strP), ...
        '�',sMYName, ...
        sprintf( ...
        'f_%s(�) - %.2f�f_res(�) = %s׵ + %s', ... %TODO check symbols
        sMYName, sqrt(pi/2)*z, ...
        num2str(polyLLoa(1)), num2str(polyLLoa(2)) ...
        ), ...
        @dcXY);
    
    % lower LOA text
    text(xl(2)-padding/2,yLLoa(2)+sign(polyLLoa(1))*padding/2, ...
            sprintf( ...
            '$f_%s(\\mu)-%.2f\\sqrt{\\pi/2}f_{res}(\\mu)$',sMYName,z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
    
    % upper LOA line
    yULoa = polyval(polyULoa,xlim);
    regLineULoa = line(xlim,yULoa);
    regLineULoa.Color = 'k';
    regLineULoa.UserData = dcStruct( ...
        sprintf('upper %s%% LOA',strP), ...
        '�', sMYName, ...
        sprintf( ...
        'f_%s(�) + %.2f�f_res(�) = %s׵ + %s', ... %TODO check symbols
        sMYName, sqrt(pi/2)*z, ...
        num2str(polyULoa(1)),num2str(polyULoa(2)) ...
        ), ...
        @dcXY);
    
    % upper LOA text
    text(xl(2)-padding/2,yULoa(2)+sign(polyULoa(1))*padding/2, ...
            sprintf( ...
            '$f_%s(\\mu)+%.2f\\sqrt{\\pi/2}f_{res}(\\mu)$',sMYName,z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
    
    % adjust y limits
    yl = ylim;
    padding = range(yl)/20; % distance to enhance visibility
    yl = [ ...
        min(yl(1), min(yLLoa)-padding), ...
        max(yl(2), max(yULoa)+padding) ...
        ];
    ylim(yl)
end

% plot least-squares line
if doPlotLS
    lsMuYLine = refline;
    lsMuYLine.Color = 'r';
    lsMuYLine.DisplayName = ...
        'least-squares';
    lsMuYLine.UserData = dcStruct( ...
        ['least-squares line of mean and ' sMYLongName], ...
        [], [], [], ... %TODO add parameters from polyfit like in correlation plot
        @dcXY);
    legEntries(end+1) = lsMuYLine;
end

% axes labels
xlabel(sprintf('mean \\it� = (M_1+M_2)/2')) %TODO check �
strYLabel = sprintf('%s \\it%s = M_1%sM_2',sMYLongName,sMYName,strYFun);
ylabel(strYLabel)
title(sprintf( ...
    ['Mean-%s plot of (%u observations):\n' ...
    ' \\rm\\itM_1\\rm: %s\n \\itM_2\\rm: %s'], ...
    sMYLongName,n,xName,yName))

% legend
legend(legEntries,'Location','SouthWest')

% reorder plot children
axMY.Children = axMY.Children([end 1:end-1]);

% % equalise axes
% axis equal

% make sure title is in visible figure area
op = axMY.OuterPosition;
axMY.OuterPosition(4) = min(op(4)+op(2),1) - op(2);
end