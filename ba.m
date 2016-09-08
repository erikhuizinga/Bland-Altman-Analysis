function varargout = ba(varargin)
%BA Bland-Altman Analysis of agreement between measurement methods.
%   Bland-Altman Analysis is a statistical method published by J. Martin
%   Bland and Douglas G. Altman in the 1980s and further developed later
%   on. It is used to determine the agreement between two measurement
%   methods that measure the same quantity. It is a popular method in
%   biostatistics and chemistry.
%
%   Syntax
%   s = BA(x,y) performs Bland-Altman Analysis on x and y, which are data
%   from two measurement methods for the same quantity respectively. x and
%   y must be numeric vectors of the same length. The calculations are done
%   at a significance level of alpha = 0.05. Output s is a structure
%   containing multiple fields with descriptive statistics about the
%   agreement of x and y. For more details on s see section Output below.
%
%   s = BA(x,y,alpha) specifies the significance level to calculate the
%   limits of agreement and confidence intervals with. alpha must be a
%   scalar in the interval [0,1]. If alpha is not specified a value of 0.05
%   is used by default to calculate 95% limits of agreement and confidence
%   intervals.
%
%   s = BA(__,Name,Value) specifies additional options using one or more
%   name-value pair arguments, in addition to any of the input arguments in
%   the previous syntaxes. For example, you can specify to create the
%   mean-difference plot using the 'PlotMeanDifference' name-value pair
%   argument.
%
%   BA(__) can be used to plot the data without returning an output
%   argument.
%
%   __ = BA(f,__) specifies the figure(s) f in which to create the plots
%   specified with the corresponding name-value pairs. The number of
%   figures in f must equal one or the number of specified plots.
%
%   __ = BA(ax,__) specifies the (array of) axes in which to create the
%   plots specified with the corresponding name-value pairs. The number of
%   axes in ax must equal the number of specified plots.
%
%   Examples
%   See and run the ba1999demo.m script for examples of the syntax of BA
%   used with data from the 1999 article by Bland and Altman. Calling BA
%   without input arguments also runs the demonstration script.
%
%   Name-Value Pair Arguments
%   Specify optional comma-separated pairs of Name,Value arguments to
%   access various options. Name is the argument name and Value is the
%   corresponding value. Name must appear inside single quotes (' '). You
%   can specify several name and value pair arguments in any order as
%   Name1,Value1,...,NameN,ValueN.
%   Example: 'XName','X','YName','Y'
%
%   'XName': Name of x variable
%   inputname of input argument x (default) | string
%   Name of x variable, specified as a string. 'XName' is used in the plot
%   titles.
%   Example: 'XName','X' sets the first measurement's name to 'X'.
%
%   'YName': Name of y variable
%   inputname of input argument y (default) | string
%   Name of y variable, specified as a string. 'YName' is used in the plot
%   titles.
%   Example: 'YName','Y' sets the second measurement's name to 'Y'.
%
%   'Exclude': Observation pairs to exclude
%   [] (default) | logical indices | numeric indices
%   Observation pairs to exclude, specified as logical or numeric indices
%   to index into x and y. The specified elements are removed from x and y
%   before any calculations or plots.
%   Example: 'Exclude',[1, 3, 4] excludes elements 1, 3 and 4 from x and y.
%   Example: 'Exclude',[0 0 1 0 1 1 0 0 1] excludes the true elements from
%   x and y.
%
%   'PlotMeanDifference': Create mean-difference plot
%   false (default) | true
%   Create the mean-difference plot if the specified value is true. The
%   mean-difference plot is a scatter plot of the difference between
%   observations versus their mean. Specifying the 'PlotAll' Name-Value
%   pair argument as true creates the mean-difference plot, regardless of
%   the 'PlotMeanDifference' value.
%
%   'PlotCorrelation': Create correlation plot
%   false (default) | true
%   Create the correlation plot if the specified value is true. The
%   correlation plot is a scatter plot of x and y. Specifying the 'PlotAll'
%   Name-Value pair argument as true creates the correlation plot,
%   regardless of the 'PlotCorrelation' value.
%
%   'PlotAll': Create all plots
%   false (default) | true
%   Create mean-difference and correlation plots if the specified value is
%   true. Setting 'PlotAll' to true overrides any value given to the
%   'PlotMeanDifference' and 'PlotCorrelation' Name-Value pair arguments.
%   However, setting it to false does not override the individual plot
%   Name-Value pair arguments.
%
%   'PlotStatistics': Add statistics to the created plots
%   'none' (default) | 'basic' | 'extended'
%   Add statistics to the created plots, specified as 'basic' or
%   'extended'. 'basic' specifies a basic set of statistics to add.
%   'extended' adds a more extended set of statistics to the plots. The
%   following statistics are added to the plots.
%   The basic set adds labelled lines for the limits of agreement to the
%   mean-difference plot. It also adds the Spearman rank correlation
%   coefficient to the legend of the mean difference plot. Furthermore, the
%   Pearson correlation coefficient is added to the legend of the
%   correlation plot.
%   The extended set adds the statistics of the basic set. Additionally,
%   the confidence intervals of the mean difference and limits of agreement
%   are plotted as error bars in the mean-difference plot. The extended set
%   does not add statistics other than the basic set to the correlation
%   plot.
%   If no plots are created, the 'PlotStatistics' value is ignored.
%
%   'PlotLeastSquares': Add a least-squares line to the created plots
%   false (default) | true
%   Add a least-squares line to the created plots if the specified value is
%   true. The least-sqares line describes the simple linear relationship
%   between x and y in the correlation plot and the mean and difference in
%   the mean-difference plot.
%
%   Output
%   The only output argument s is optional. It is a scalar structure
%   containing multiple fields with descriptive statistics about the
%   agreement of x and y. s contains the following fields:
%   
%   muD: the mean difference between x and y, also called the bias.
% 
%   muDCI: the 95% (default, depending on alpha) confidence interval of
%   the mean difference.
% 
%   loa: the 95% (default, depending on alpha) limits of agreement, a 2
%   element vector. The first element is the lower limit of agreement, the
%   second is the upper.
% 
%   loaCI: the 95% (default, depending on alpha) confidence interval of the
%   limits of agreement, a 2x2 matrix. The first column corresponds to
%   lower limit of agreement, the second to the upper limit. The first and
%   second row correspond to the lower and upper confidence interval bound
%   respectively.
% 
%   sD: the standard deviation of the differences.
% 
%   rSMuD: the Spearman rank correlation between mean and difference.
% 
%   pRSMuD: the p-value of the Spearman rank correlation for testing the
%   hypothesis of no correlation against the alternative that there is a
%   nonzero correlation. %TODO also output rhoXY,pRhoXY
%
%   About
%   This MATLAB function is an implementation of the methods in the 1999
%   article by Bland and Altman:
%   Bland & Altman 1999 - Measuring agreement methods in comparison studies
%   http://smm.sagepub.com/content/8/2/135.abstract
%
%   You might not have access to this article. Access it through your
%   institution's library or buy it.
%
%   The article comprises of 5 methodological sections (sections 2-6). The
%   current version of this MATLAB file (8 september 2016) implements the
%   first methodological section. More sections will be added in the
%   future.
%
%   The demonstration script `ba1999demo.m` is an implementation of the
%   calculations done by Bland and Altman in the article. Their article
%   contains a number of example data sets, which they use in their
%   methods. The demonstration script illustrates the same results and the
%   syntax used to obtain them.
%
%   See also BA1999DEMO

%   Copyright (C) 2016 Erik Huizinga, huizinga.erik@gmail.com
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% inputs
% demo if no input arguments
if ~nargin, ba1999demo, return, end

in = varargin;

% check if first input is (an array of) handles
if all(isgraphics(in{1}))
    h = in{1};
    in = in(2:end);
    ixname = 2;
    iyname = 3;
else
    h = [];
    ixname = 1;
    iyname = 2;
end
% default axes variables
axMD = [];
axC = [];

% prepare parser
p = inputParser;
p.addRequired('x',@isnumeric)
p.addRequired('y',@isnumeric)
p.addOptional('a',.05,@isnumeric)
p.addParameter('XName',inputname(ixname),@ischar)
p.addParameter('YName',inputname(iyname),@ischar)
p.addParameter('PlotAll',false,@validatelogical)
p.addParameter('PlotMeanDifference',false,@validatelogical)
p.addParameter('PlotCorrelation',false,@validatelogical)
p.addParameter('Exclude',[],@validatelogical)
p.addParameter('PlotStatistics','none',@ischar)
p.addParameter('PlotLeastSquares',false,@validatelogical)

% parse inputs
parse(p,in{:})
s2v(p.Results);

%% validate and preprocess inputs
% x and y: measurements of two methods
if numel(x)~=numel(y)
    error 'Number of elements in x and y must be equal.'
end
% force column vectors
if isvector(x), x = x(:); end
if isvector(y), y = y(:); end

% alpha: significance level
if ~isscalar(a) && a<0 && a>1
    error 'alpha must be a scalar in the interval [0,1].'
end

% xName and yName
if ~iscellstr(XName), XName = cellstr(XName); end
if ~iscellstr(YName), YName = cellstr(YName); end
xName = strjoin(XName,', ');
yName = strjoin(YName,', ');

% doAllPlots
doPlotAll = logical(PlotAll);
if doPlotAll
    doPlotMD = true; % do mean-difference plot
    doPlotC = true; % do correlation plot
else
    doPlotMD = logical(PlotMeanDifference);
    doPlotC = logical(PlotCorrelation);
end

% validate number and type of handles in h for the requested plots
doPlot = [doPlotMD doPlotC];
if any(doPlot)
    nPlot = nnz(doPlot);
    if isempty(h)
        figure
        ax = gobjects(0);
    else
        nh = numel(h);
        if isaxes(h)
            if nh==nPlot
                ax = h;
            else
                error(['The number of requested plots (' num2str(nPlot) ...
                    ') is not equal to the number of axes (' ...
                    num2str(nh) ') in the first input argument.'])
            end
        else % h must be a figure or contain nplot figures
            if nh==1
                figure(h) % error if h is not a figure
                ax = gobjects(0);
            elseif nh==nPlot
                for n = nh:-1:1
                    figure(h(n))
                    ax(n) = axes;
                end
            elseif all(isfigure(h))
                error(['The number of requested plots (' num2str(nPlot) ...
                    ') is not equal to the number of figures (' ...
                    num2str(nh) ') in the first input argument.'])
            else
                error(['The first optional input argument must be an ' ...
                    '(array of) axes or figure handle(s).'])
            end
        end
    end
    
    % create axes if not yet existent
    if isempty(ax)
        if nPlot>1
            for n = nPlot:-1:1
                ax(n) = subplot(1,nPlot,n);
            end
        else
            ax = axes;
        end
    end
    
    % store axes in plot order
    if doPlotMD, axMD = ax(1); ax(1) = []; end
    if doPlotC, axC = ax(1); end % ax(1) = []; end
end

% exclude samples
lex = false(size(x));
lex(Exclude) = true;
x(lex) = [];
y(lex) = [];

% doPlotStats
% doPlotStats = logical(PlotStatistics);

% doPlotAddStats
% doPlotAddStats = logical(PlotAdditionalStatistics);
switch lower(PlotStatistics)
    case 'none'
        doPlotBasicStats = false;
        doPlotExtStats = false;
    case 'basic'
        doPlotBasicStats = true;
        doPlotExtStats = false;
    case 'extended'
        doPlotBasicStats = true;
        doPlotExtStats = true;
    otherwise
        error 'Unknown value for the ''PlotStatistics'' name-value pair.'
end

% doPlotLS
doPlotLS = logical(PlotLeastSquares);

% % dobaa4r (do BAA for repeated measurements)
% doBAA4r = logical(Repeated);
% if doBAA4r && isvector(x)
%     warning(['x and y must be matrices with every row an observation ' ...
%         'and every column the observed data. Doing normal ' ...
%         'Bland-Altman analysis instead.'])
%     doBAA4r = false;
% end

%% Bland-Altman analysis
% if doBAA4r, fba = @baa4r; else fba = @baa; end
% out = fba(x, xName, y, yName, a, ...
%     doPlotMD, axMD, ...
%     doPlotC, axC, ...
%     doPlotStats, doPlotAddStats, doPlotLS);
out = baa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtStats, doPlotLS);

%% output
if nargout, varargout = out; end
end

function s2v(s)
%#ok<*NODEF>
fn = fieldnames(s);
for f = fn(:).'
    assignin('caller',f{:},s.(f{:}))
end
end

function tf = validatelogical(L)
%VALIDATELOGICAL True for arrays that can be converted to logical.
%   VALIDATELOGICAL(L) returns logical 1 (true) if L is a logical or
%   numeric array and logical 0 (false) otherwise.
%
%   See also LOGICAL.

try tf = islogical(logical(L)); catch, tf = false; end
end

function varargout = baa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtStats, doPlotLS)
%% preparation
% multiFig = false;
if doPlotMD || doPlotC
    % if doMDPlot, f(1) = axMD.Parent; end
    % if doCPlot, f(2) = axC.Parent; end
    ax = [axMD;axC];
    f = get(ax,'Parent');
    if iscell(f), f = vertcat(f{:}); end
    f = unique(f);
    % f can be the handle to one or more figures
    % if numel(f)>1, multiFig = true; end
else
    f = [];
end

% keep only values that can be used in calculations
lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
n = nnz(lok);
xok = x(lok);
yok = y(lok);

%% calculation
% mean statistics
muXY = mean([xok,yok],2);

% difference statistics
d = xok-yok; % difference
muD = mean(d); % mean of difference

sD = std(d); % s_d in article
% varD = var(d); % s_d^2 in article

varMuD = sD^2/n; % variance of muD (SE^2) (p. 141)
seMuD = sqrt(varMuD); % standard error of muD (p. 142)

% varSD = sD^2/2/(n-1); % approximated variance of sD (article p. 141)
% Instead, use exact formula for unbiased estimated variance
% Source: http://stats.stackexchange.com/a/28567/80486
% sSD = sD * gamma((n-1)/2) / gamma(n/2) * ...
%     sqrt( (n-1)/2 - ( gamma(n/2) / gamma((n-1)/2) )^2 ); 
gammafrac = gamma((n-1)/2) / gamma(n/2);
% if ~isfinite(gammafrac) % true for large n
    % % approximate using gamma(a+b)/gamma(a) ~ a^b
    % % Source: https://en.wikipedia.org/w/index.php?title=Gamma_function&oldid=737220343#General
    % % compare:
    % % figure
    % % n = 1:500;
    % % g1 = gamma((n-1)/2)./gamma(n/2);
    % % g2 = ((n-1)/2).^(-1/2);
    % % g3 = sqrt(2./(n-1));
    % % plot(n,g1,n,g2,n,g3,n,g1-g2,n,g1-g3,n,g2-g3)
    % % legend g1 g2 g3 g1-g2 g1-g3 g2-g3
    % gammafrac = sqrt(2/(n-1)); % same as: gammafrac = ((n-1)/2).^(-1/2);
% end
sSD = sD * gammafrac * sqrt( (n-1)/2 - gammafrac^-2 );
varSD = sSD^2; % unbiased estimate of variance of s_d

% limits of agreement statistics
p = 1-a/2;
z = Ninv(p); % inverse normal distribution at p = 1-alpha/2
loa = muD + z*sD*[-1 1]; % limits of agreement (LOA)

% confidence intervals (CI) for muD and loa
t = Tinv(p,n-1);
varLoa = varMuD + z^2*varSD; % article: bottom of p. 141
seLoa = sqrt(varLoa); % standard error of the LOA
eLoa = t*seLoa; % LOA error
eMuD = t*seMuD; % muD error
muDCI = muD + eMuD*[-1 1];
loaCI = [loa;loa] + eLoa*[-1 -1;1 1];
% loaCI is a 2x2 matrix. Every column contains the CI for a LOA: the first
% column corresponds to loa(1), the second to loa(2). The first and second
% row correspond to the lower and upper CI bound respectively.
% LOA = [loaCI(1,:);loa;loaCI(2,:)]; % optional 3x2 matrix form

% mean difference correlation statistics
[rSMuD,pRSMuD] = corr(muXY,d,'type','Spearman'); %TODO make independent of stats toolbox?

% linear regression and correlation %TODO add polyfit for muXY and d
[polyXY,statsPXY] = polyfit(xok,yok,1);
[rhoXY,pRhoXY] = corrcoef(xok,yok);
rhoXY = rhoXY(1,2);
pRhoXY = pRhoXY(1,2);

% linear regression statistics
resPXY = yok - polyval(polyXY,xok); % residual
ssePXY = sum(resPXY.^2); % or rss/ssr
msePXY = ssePXY/statsPXY.df;
% compare with the same calculation (CF toolbox required):
% [~,gof] = fit(xok,yok,'Poly1');
% gof.sse-sse % equal within double precision
% gof.rmse-sqrt(mse) % idem
% gof.dfe-s.df % equal

%% graphics
% correlation plot
if doPlotC
    % preparation
    axes(axC)
    legEntries = gobjects(0);
    
    % plot y against x
    sC = scatter(xok,yok);
    sC.UserData = dcStruct([],'M1','M2',[],@dcXY);
    axis tight
    axis equal
    
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
        XYLine = refline;
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
    eqLine = refline(1,0); % 45Â° degree, because of axis equal
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
end

% mean-difference plot
if doPlotMD
    % preparation
    axes(axMD)
    legEntries = gobjects(0);
    
    % mean-difference plot
    sMD = scatter(muXY,d);
    sMD.UserData = dcStruct([],'µ','d',[],@dcXY); %TODO check µ
    axis equal
    xl = xlim;
    
    % plot statistics
    if doPlotBasicStats
        % prepare axes
        hold(axMD,'on')
        padding = range(ylim)/20; % distance to enhance visibility
        yl = ylim;
        yl = [ ...
            min(yl(1), loaCI(1,1)) - padding, ...
            max(yl(2), loaCI(2,2)) + padding ...
            ];
        ylim(yl)
        
        % add correlation to legend
        if pRSMuD<1e-4
            % as recommended by BMJ 1996;312:572
            % http://dx.doi.org/10.1136/bmj.312.7030.572
            strPRSMuD = sprintf('\\itp\\rm < 0.0001');
        else
            strPRSMuD = sprintf('\\itp\\rm = %.2f',pRSMuD);
        end
        sMD.DisplayName = sprintf( ...
            '\\itr_s\\rm = %.2f (%s)',rSMuD,strPRSMuD);
        legEntries(end+1) = sMD;
        
        % lower LOA line
        lLine = refline(0,loa(1));
        lLine.Color = 'k';
        lLine.UserData = dcStruct([],[],'lower LOA', ...
            [num2str(100*(1-a)) '% CI = [' num2str(loaCI(1,1)) ', ' num2str(loaCI(2,1)) ']'], ...
            @dcXY);
        text(xl(2)-padding/2,loa(1)+padding/2,sprintf( ...
            '$\\overline{d}-%.2fs_d$',z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % upper limit of agreement line
        uLine = refline(0,loa(2));
        uLine.Color = 'k';
        uLine.UserData = dcStruct([],[],'upper LOA', ...
            [num2str(100*(1-a)) '% CI = [' ...
            num2str(loaCI(1,2)) ', ' num2str(loaCI(2,2)) ']'], ...
            @dcXY);
        text(xl(2)-padding/2,loa(2)+padding/2, ...
            sprintf( ...
            '$\\overline{d}+%.2fs_d$',z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % mean difference line
        muDLine = refline(0,muD); % mean difference
        muDLine.Color = 'k';
        text(xl(2)-padding/2,muD+padding/2,'$\overline{d}$', ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        muDLine.UserData = dcStruct([],[],'mean difference', ...
            [num2str(100*(1-a)) '% CI = [' ...
            num2str(muDCI(1)) ', ' num2str(muDCI(2)) ']'], ...
            @dcXY);
        
        % plot additional statistics
        if doPlotExtStats
            % lower LOA errorbar
            lE = errorbar(xl(2),loa(1),eLoa,'r');
            lE.UserData = lLine.UserData;
            
            % mean difference errorbar
            muDE = errorbar(xl(2),muD,eMuD,'r');
            muDE.UserData = muDLine.UserData;
            
            % upper LOA errorbar
            uE = errorbar(xl(2),loa(2),eLoa,'r');
            uE.UserData = uLine.UserData;
        end
    end
    
    % plot least-squares line
    if doPlotLS
        % least-squares line
        lsMuDLine = refline;
        lsMuDLine.Color = 'r';
        lsMuDLine.DisplayName = ...
            'least-squares';
        lsMuDLine.UserData = dcStruct( ...
            'least-squares line of mean and difference', ...
            [], [], [], ... %TODO add parameters from polyfit like in correlation plot
            @dcXY);
        legEntries(end+1) = lsMuDLine;
    end
    
    % zero difference line
    line0 = refline(0,0);
    line0.Color = [.75 .75 .75];
    line0.LineStyle = '--';
    line0.DisplayName = 'line of equality';
    line0.UserData = dcStruct(line0.DisplayName,[],[],[],@dcXY);
    legEntries(end+1) = line0;
    
    % axes labels
    xlabel(sprintf('mean \\itµ = (M_1+M_2)/2')) %TODO check µ
    ylabel(sprintf('difference \\itd = M_1-M_2'))
    title(sprintf(['Mean-difference plot of (%u observation pairs):\n' ...
        ' \\rm\\itM_1\\rm: %s\n \\itM_2\\rm: %s'],n,xName,yName))
    
    % legend
    legend(legEntries,'Location','SouthWest')
    
    % reorder plot children
    axMD.Children = axMD.Children([end 1:end-1]);
end

%% set data cursor update function for figure(s)
for f = f(:).'
    dc = datacursormode(f);
    dc.UpdateFcn = @dcUpdateFcn;
    dc.SnapToDataVertex = 'off';
    dc.Enable = 'on';
end

%% output
out.muD = muD;
out.muDCI = muDCI;
out.loa = loa;
out.loaCI = loaCI;
out.sD = sD;
out.rSMuD = rSMuD;
out.pRSMuD = pRSMuD;
varargout = {{out}};
end

function x = Ninv(p)
% inverse normal
x = -sqrt(2).*erfcinv(2*p);
end

function x = Tinv(p,n)
% inverse t
assert(p>=.75 & p<=1,'p must satisfy 0.75<=p<=1.')
if n<=1
    error 'Not implemented for n<=1.'
elseif n<1000
    bii = betaincinv(2*abs(p-1/2), n/2, 1/2, 'upper');
    x = sqrt(n*(1-bii)/bii);
else % n>=1000
    % Reference:
    % Abramowitz & Stegun, 10th ed. 1972, formula 26.7.5
    xp = Ninv(p);
    g1 = @(x) (x^3+x)/4;
    g2 = @(x) (5*x^5+16*x^3+3*x)/96;
    g3 = @(x) (3*x^7+19*x^5+17*x^3-15*x)/384;
    g4 = @(x) (79*x^9+776*x^7+1482*x^5-1920*x^3-945*x)/92160;
    x = xp + g1(xp)/n + g2(xp)/n^2 + g3(xp)/n^3 + g4(xp)/n^4;
end
end

function tf = isfigure(f)
tf = strcmp({f.Type},'figure');
end

function tf = isaxes(ax)
tf = strcmp({ax.Type},'axes');
end

function str = dcUpdateFcn(~,e)
x = e.Position(1);
y = e.Position(2);
data = e.Target.UserData;
str = data.dcFun(data,x,y);
if ~isempty(data.prefix), str = [data.prefix, str]; end
if ~isempty(data.suffix), str = [str, data.suffix]; end
end

function s = dcStruct(prefix,xName,yName,suffix,dcFun)
s.prefix = prefix;
s.xName = xName;
s.yName = yName;
s.suffix = suffix;
s.dcFun = dcFun;
end

function str = dcXY(data,x,y)
str = {};
if ~isempty(data.xName), str = [str, [data.xName ': ' num2str(x)]]; end
if ~isempty(data.yName), str = [str, [data.yName ': ' num2str(y)]]; end
end