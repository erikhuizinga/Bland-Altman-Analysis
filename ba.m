function varargout = ba(varargin)
%BA Bland-Altman Analysis of agreement between measurement methods.
%   Bland-Altman Analysis is a statistical method published by J. Martin
%   Bland and Douglas G. Altman in the 1980s and further developed later
%   on. It is used to determine the agreement between two measurement
%   methods that measure the same quantity. It is a popular method in
%   biostatistics and chemistry.
%
%   Syntax
%   stats = BA(x, y) performs Bland-Altman Analysis (BAA) on x and y, which
%   are data from two measurement methods for the same quantity
%   respectively. x and y can be of various classes and shapes:
%    - If x and y are vectors, regular BAA is performed. Every element in x
%      corresponds to the element in y at the same index. These pairs are
%      individual observations on individual subjects.
%    - If x and y are matrices, BAA for repeated measurements is performed.
%      This means multiple measurements have been acquired per subject. The
%      number of elements in x and y must be equal. Every row of x and y
%      corresponds to the subjects, every column to the repeated
%      observations.
%   For more information about BAA for repeated measurements, see section
%   Bland-Altman Analysis for repeated measurements below. The calculations
%   are done at a significance level of alpha = 0.05. Output stats is a
%   structure containing multiple fields with descriptive statistics about
%   the agreement of x and y. For more details on the fields in stats see
%   section Output below.
%
%   stats = BA(x, y, alpha) specifies the significance level to calculate
%   the limits of agreement and confidence intervals with. alpha must be a
%   scalar in the interval [0, 1]. If alpha is not specified a value of
%   0.05 is used by default to calculate 95% limits of agreement and
%   confidence intervals.
%
%   stats = BA(__, Name, Value) specifies additional options using one or
%   more Name-Value pair arguments, in addition to any of the input
%   arguments in the previous syntaxes. For example, you can specify to
%   create the mean-difference plot using the 'PlotMeanDifference'
%   Name-Value pair argument.
%
%   BA(__) can be used to create graphs without returning an output
%   argument.
%
%   __ = BA(f, __) specifies the figure(s) f in which to create the graphs
%   specified with the corresponding Name-Value pairs. The number of
%   figures in f must equal one or the number of specified plots. If f is
%   the handle to only one figure, a number of subplots is created in the
%   figure equal to the number of graphs requested.
%
%   __ = BA(ax, __) specifies the (array of) axes in which to create the
%   requested graphs with the corresponding Name-Value pairs. The number of
%   axes in ax must equal the number of graphs requested.
%
%   Examples
%   See and run the bademo.m script for examples of the syntax of BA used
%   with data from the 1999 article by Bland and Altman. Calling BA without
%   input arguments also runs the demonstration script.
%
%   Name-Value Pair Arguments
%   Specify optional comma-separated pairs of Name,Value arguments to
%   access various options. Name is the argument name and Value is the
%   corresponding value. Name must appear inside single quotes (' '). You
%   can specify several name and value pair arguments in any order as
%   Name1, Value1, ..., NameN, ValueN.
%   Example: 'XName', 'X', 'YName', 'Y'
%
%   'XName': Name of x variable
%   inputname() of input argument x (default) | string
%   Name of x variable, specified as a string. 'XName' is used in the plot
%   titles.
%   Example: 'XName', 'X' sets the first measurement's name to 'X'.
%
%   'YName': Name of y variable
%   inputname() of input argument y (default) | string
%   Name of y variable, specified as a string. 'YName' is used in the plot
%   titles.
%   Example: 'YName', 'Y' sets the second measurement's name to 'Y'.
%
%   'Exclude': Subjects to exclude
%   [] (default) | logical indices | numeric indices
%   Subjects to exclude, specified as logical or numeric indices into x and
%   y. The specified rows are removed from x and y before any calculations
%   are done or graphs are created.
%   Example: 'Exclude', [1, 3, 4] excludes rows 1, 3 and 4 from x and y.
%   Example: 'Exclude', [0 0 1 0 1 1 0 0 1] excludes the true rows from
%   x and y. Note the logical vector needs not be a column vector.
%
%   'Transform': Function to transform data with
%   @(x) x (default) | function handle
%   Function to transform data with before further analysis, specified as a
%   function handle of one variable. By default, no transformation is
%   performed. The function handle should accept a vector input. Bland and
%   Altman suggest in their 1999 article (see p. 144) only the logarithmic
%   transformation should be used, i.e. specify 'Transform', @log. Other
%   transforms are not easily relatable to the actual measurements, hence
%   their recommendation. Note that no backtransformation is performed on
%   the statistics in the optional output argument, i.e. they are the
%   transformed statistics.
%   Example: 'Transform', @log transforms x to log(x) and y to log(y).
%
%   'PlotMeanDifference': Create mean-difference graph
%   false (default) | true
%   Create the mean-difference graph if the specified value is true. The
%   mean-difference graph is a scatter plot of the difference between
%   observations versus their mean. Specifying the 'PlotDefault' Name-Value
%   pair argument as true creates the mean-difference plot, regardless of
%   the 'PlotMeanDifference' value.
%
%   'PlotMeanRatio': Create mean-ratio graph
%   false (default) | true
%   Create the mean-ratio graph if the specified value is true. The
%   mean-ratio graph is a scatter plot of the ratio between observations
%   versus their mean. If the mean-ratio graph is created and an output
%   arguments is specified, the output argument contains a field called
%   ratio with the ratio statistics.
%   
%   'PlotMeanSD': Create mean-standard deviation graph
%   false (default) | true | string
%   Create the mean-standard deviation graph if the specified value is not
%   false or 'none'. The mean-standard deviation graph is a scatter plot of
%   a standard deviation versus the mean. The plotted standard deviation
%   depends on the specified value. The options for the value are listed
%   below. The standard deviation of the difference, ratio or of x and y
%   can be plotted. If the standard deviations of x and y are plotted, two
%   axes are required. The following values are allowed:
%   true | 'difference' | 'single' | 'joint': plot the standard deviation
%   of the difference versus the mean.
%   'ratio': plot the standard deviation of the ratio versus the mean.
%   'input' | 'both' | 'separate': plot the standard deviation of x and y
%   in separate axes versus the mean.
%   false (default) | 'none': no mean-standard deviation graph is created.
%
%   'PlotCorrelation': Create correlation graph
%   false (default) | true
%   Create the correlation graph if the specified value is true. The
%   correlation graph is a scatter plot of x and y. Specifying the
%   'PlotAll' Name-Value pair argument as true creates the correlation
%   plot, regardless of the 'PlotCorrelation' value. If the correlation
%   graph is created and an output arguments is specified, the output
%   argument contains a field called correlation with the correlation
%   statistics.
%
%   'PlotDefault': Create default graphs
%   false (default) | true
%   Create mean-difference and correlation graphs if the specified value is
%   true. Setting 'PlotDefault' to true overrides any value given to the
%   'PlotMeanDifference' and 'PlotCorrelation' Name-Value pair arguments.
%   However, setting 'PlotDefault' to false does not override the
%   individual Name-Value pair arguments. Setting 'PlotDefault' to true and
%   specifying the optional output argument adds the correlation field to
%   the output argument.
%
%   'PlotStatistics': Add statistics to the created plots
%   'none' (default) | 'basic' | 'extended'
%   Add statistics to the created plots, specified as 'none', 'basic',
%   'extended' or 'regression'. 'none' specifies no statistics to be added
%   to the graphs. 'basic' specifies a basic set of statistics to add.
%   'extended' adds a more extended set of statistics. 'regression' adds
%   regression lines to the graphs.The following statistics are added to
%   the plots. The basic set adds labelled lines for the limits of
%   agreement to the mean-statistic graphs. It also adds the Spearman rank
%   correlation coefficient to the legend of these graphs. Furthermore, the
%   Pearson correlation coefficient is added to the legend of the
%   correlation plot. The extended set adds the statistics of the basic
%   set. Additionally, the confidence intervals of the bias and limits of
%   agreement are plotted as error bars in the mean-statistic graphs. The
%   extended set does not add statistics other than the basic set to the
%   correlation plot. The regression statistics comprise of the extended
%   set, except that the bias and limits of agreement lines are no longer
%   constant with respect to the mean. The variable on the vertical axis is
%   regressed on the mean, resulting in the possibility of non-constant
%   lines. If no plots are created, the 'PlotStatistics' value is ignored.
%
%   'ConstantResidualVariance': Assume constant residual variance
%   false (default) | true
%   Assume constant residual variance in the simple linear regression
%   performed if the 'PlotStatistics', 'regression' Name-Value pair
%   argument is specified. This means the upper and lower limits of
%   agreement lines will have the same slope as the bias line. This
%   assumption holds if the slope of the upper and lower limits of
%   agreement do not differ significantly from the slope of the bias
%   regression line.
%
%   'ConstantTrueValue': Assume the true value is constant
%   true (default) | false
%   Assume the true value being measured by x and y is constant. This
%   assumption is true for repeated measurements of the same constant
%   quantity. For example, repeated measurements of a subject's blood
%   pressure within the same minute could be assumed to be constant. If
%   this is to be assumed in the calculations, specify the 
%   'ConstantTrueValue', true pair argument. If the repeated measurements
%   are of a quantity that changes from measurement to measurement, this
%   assumption does not hold. For example, blood pressure measured on
%   multiple occasions throughout a day is expected to change over time,
%   hence over repeated measurements as well. If a varying true value is
%   assumed, specify the 'ConstantTrueValue', false pair argument.
%
%   Output
%   The only output argument, stats, is optional. It is a scalar structure
%   containing multiple fields with descriptive statistics about the
%   agreement of x and y. The number of fields in stats varies depending on
%   the input arguments. By default stats contains three fields, being
%   difference, x and y. Additional fields in stats are ratio and
%   correlation. They exist depending on the requested graphs. If the
%   mean-ratio graph is requested using the 'PlotMeanRatio' Name-Value pair
%   argument, then stats is returned with the ratio field. If the
%   correlation plot is requested using the 'PlotCorrelation' or
%   'PlotDefault' Name-Value pair arguments, then stats is returned with
%   the correlation field. Note that the difference field is returned
%   regardless of the 'PlotMeanDifference' Name-Value pair argument.
%
%   Fields difference and ratio are described together below, after which
%   fields x and y are described together too. Lastly, field correlation is
%   described.
%
%   The difference and ratio fields are structures themselves, containing
%   the statistics about the differences and ratios repectively between
%   observations in x and y. The fields in difference and ratio are
%   described below (the word statistic refers to either difference or
%   ratio):
%
%   mu: the mean statistic between x and y, also called the bias.
%   Example: stats.difference.mu is the mean difference.
%
%   muCI: the 95% (default, depending on alpha) confidence interval of
%   the mean statistic.
%   Example: stats.ratio.muCI is the confidence interval of the mean ratio.
%
%   loa: the 95% (default, depending on alpha) limits of agreement, a 2
%   element vector. The first element is the lower limit of agreement, the
%   second is the upper.
%   Example: stats.difference.loa(1) is the lower limit of agreement of the
%   differences.
%
%   loaCI: the 95% (default, depending on alpha) confidence interval of the
%   limits of agreement, a 2x2 matrix. The first column corresponds to
%   lower limit of agreement, the second to the upper limit. The first and
%   second row correspond to the lower and upper confidence interval bound
%   respectively.
%   Ecample: stats.ratio.loaCI(:, 2) is the confidence interval of the
%   upper limit of agreement of the ratios.
%
%   s: the standard deviation of the statistic.
%   Example: stats.difference.s is the standard deviation of the
%   differences.
%
%   rSMu: the Spearman rank correlation between mean and statistic.
%   Example: stats.ratio.rSMu is the Spearman rank correlation between mean
%   and the ratios.
%
%   pRSMu: the p-value of the Spearman rank correlation for testing the
%   hypothesis of no correlation against the alternative that there is a
%   nonzero correlation.
%   Example: stats.difference.pRSMu is the p-value of the Spearman rank
%   correlation between mean and difference, to test the hypothesis of no
%   correlation against the alternative of nonzero correlation.
%
%   polyMu: the polynomial coefficients of the simple linear regression of
%   the statistic on the mean. The first element of polyMu is the slope,
%   the second the intercept.
%   Example: stats.ratio.polyMu(2) is the intercept of the simple linear
%   regression line of the ratio on the mean.
%
%   msePolyMu: the mean squared error (MSE) of the simple linear regression
%   of the statistic on the mean.
%   Example: stats.difference.msePolyMu is the MSE of the simple linear
%   regression of the difference on the mean.
%
%   The x and y fields of structure stats contain statistics about inputs x
%   and y that are calculated in BAA. The fields in x and y are structures
%   themselves and have the following fields (the word input referring to x
%   or y):
%
%   varWithin: the within-subject variance of the input. varWithin equals
%   zero for BAA without repeated measurements, but can be non-zero in BAA
%   for repeated measurements.
%   Example: stats.x.varWithin is the within-subject variance of input x.
%
%   If the correlation field exists in stats, it contains the following
%   fields:
%
%   rho: the Pearson correlation coefficient between inputs x and y.
%
%   p: the p-value of the Pearson correlation coefficient between inputs x
%   and y.
%
%   poly: the polynomial coefficients of the simple linear regression of y
%   on x. The first element in poly is the slope, the second the intercept.
%   Example: stats.correlation.poly(1) is the slope of the simple linear
%   regression line of y on x.
%
%   polyMSE: the mean squared error (MSE) of the simple linear regression
%   of y on x.
%
%   Bland-Altman Analysis for repeated measurements
%   When performing regular limits of agreement (LOA) estimation through
%   Bland-Altman Analysis (BAA) a number of implicit assumptions exists.
%   One of these assumptions is that the individual observation pairs are
%   independent. This is the case when every observed pair comes from a
%   different subject. However, when multiple measurements are taken on the
%   same individual, this assumption migtht not hold. The measurements on
%   the same individual might depend on each other. For example, a
%   subject's blood pressure is measured every minute using a
%   sphygmomanometer. The observations of every minute will be quite
%   correlated with the previous minute, i.e. the cross-covariance is high
%   at the first (few) lag(s). Another more general example: consider 10
%   observations of a physical quantity. If these 10 observations are
%   obtained from 10 different subjects (assuming these subjects are a
%   representative sample of the population), the assumption of
%   uncorrelatedness in regular BAA will have been met. However, if these
%   10 observations would come from one subject only (as an extreme
%   example), the calculations of (mean) difference (or ratio) will depend
%   strongly on the within-subject variance of this particular subject.
%   This is less informative about the performance of the two measurement
%   methods in the population of interest and how they compare in general.
%   Instead, this analysis tells us about the performance of the methods on
%   this particular subject, which usually is not of interest in method
%   comparison studies.
%
%   About & References
%   This MATLAB function is an implementation of the methods in the 1983,
%   1986, 1999 and 2007 articles by Bland and Altman, especially these two:
%
%   [BA1999] Bland, J.M., and Altman, D.G. 1999. Measuring agreement in
%   method comparison studies. Statistical Methods in Medical Research 8:
%   135-160. doi:10.1191/096228099673819272.
%   [BA2007] Bland, J.M., and Altman, D.G. 2007. Agreement Between Methods
%   of Measurement with Multiple Observations Per Individual. Journal of
%   Biopharmaceutical Statistics 17: 571-582.
%   doi:10.1080/10543400701329422.
%
%   You might not have access to these articles. Access them through your
%   institution's library or buy access to read them.
%
%   The demonstration script bademo.m is an implementation of the
%   calculations done by Bland and Altman in their articles. Their articles
%   contain a number of example data sets, which they use in their methods
%   to demonstrate them. The demonstration script illustrates the same
%   results and the syntax used to obtain them.
%
%   See also BADEMO


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


%% Parse inputs
% Run bademo if no input arguments are given
if ~nargin, bademo, return, end


% Check if first input is (an array of) handles
in = varargin;
if all(isgraphics(in{1}))
    h = in{1};
    in = in(2 : end);
    ixname = 2;
    iyname = 3;
else
    h = [];
    ixname = 1;
    iyname = 2;
end


% Prepare parser
p = inputParser;
p.addRequired('x', @validateXY)
p.addRequired('y', @validateXY)
p.addOptional('a', .05, @isnumeric)
p.addParameter('Exclude', [], @validatelogical)
p.addParameter('XName', inputname(ixname), @ischar)
p.addParameter('YName', inputname(iyname), @ischar)
p.addParameter('PlotDefault', false, @validatelogical)
p.addParameter('PlotMeanDifference', false, @validatelogical)
p.addParameter('PlotMeanRatio', false, @validatelogical)
p.addParameter('PlotMeanSD', false, @(s) validatelogical(s)|ischar(s))
p.addParameter('PlotCorrelation', false, @validatelogical)
p.addParameter('PlotStatistics', 'none', @ischar)
p.addParameter('ConstantResidualVariance', false, @validatelogical)
p.addParameter('ConstantTrueValue', true, @validatelogical)
p.addParameter('Transform', @(x) x, ...
               @(f) isa(f,'function_handle') | ischar(f))

% Parse inputs
parse(p, in{:})
s2v(p.Results);  %#ok<*NODEF>


%% Validate and preprocess inputs
% Validate x and y: measurements of two methods
% Validate and reshape x and y for further analysis, check for repeated
% measurements
[xok, yok, doRepeated] = parseXY(x, y);


% Validate alpha: significance level
if ~isscalar(a) && a<0 && a>1
    error 'alpha must be a scalar in the interval [0,1].'
end


% Set xName and yName
if ~iscellstr(XName), XName = cellstr(XName); end
if ~iscellstr(YName), YName = cellstr(YName); end
xName = strjoin(XName, ', ');
yName = strjoin(YName, ', ');


% Validate plot arguments
[doPlotMD, axMD, doPlotMR, axMR, ...
    MSDType, doPlotMSD1, axMSD1, doPlotMSD2, axMSD2, doPlotC, axC] ...
    = validatePlotArgs(PlotDefault, PlotMeanDifference, PlotMeanRatio, ...
                       PlotMeanSD, PlotCorrelation, h);

                   
% Exclude unwanted samples
lex = false(size(xok, 1), 1);
lex(Exclude) = true;
xok(lex, :) = [];
yok(lex, :) = [];


% Set statistics set to plot
[doPlotBasicStats, doPlotExtendedStats, ...
    doPlotRegStats, doConstantRegression] ...
    = parseStatArgs(PlotStatistics, ConstantResidualVariance);


% Set transformation function
transFun = Transform;
if ischar(transFun), transFun = str2func(transFun); end
switch lower(char(transFun))
    % Detect supported transformations
    case 'log'
        if any(strcmpi(p.UsingDefaults, 'XName'))
            xName = ['Log ', xName];
        end
        if any(strcmpi(p.UsingDefaults, 'YName'))
            yName = ['Log ', yName];
        end
        
        % Transform
        xok = transFun(xok);
        yok = transFun(yok);
    otherwise
        % Do not transform
end


% Set constant true value (CTV) assumption
assumeCTV = logical(ConstantTrueValue);


%% Perform Bland-Altman analysis
out = baa( ...  % baa: Bland-Altman Analysis
           xok, xName, yok, yName, a, ...
           doPlotMD, axMD, ...
           doPlotMR, axMR, ...
           MSDType, doPlotMSD1, axMSD1, doPlotMSD2, axMSD2, ...
           doPlotC, axC, ...
           doPlotBasicStats, doPlotExtendedStats, ...
           doPlotRegStats, doConstantRegression, ...
           doRepeated, assumeCTV);


%% Set output
if nargout, varargout = out; end
end