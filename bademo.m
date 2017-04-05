%BADEMO Demonstration of Bland-Altman-Analysis, a MATLAB implementation by
%   Erik Huizinga based on the methods by Bland and Altman 1983, 1986, 1999
%   and 2007. The references below are used throughout this demo script.
%
% References in the comments of the code
% [BA1999] Bland, J.M., and Altman, D.G. 1999. Measuring agreement in method comparison studies. Statistical Methods in Medical Research 8: 135-160. doi:10.1191/096228099673819272.
% [BA2007] Bland, J.M., and Altman, D.G. 2007. Agreement Between Methods of Measurement with Multiple Observations Per Individual. Journal of Biopharmaceutical Statistics 17: 571-582. doi:10.1080/10543400701329422.
% [BSOCA1993] Bowling, L.S., Sageman, W.S., O’Connor, S.M., Cole, R., and Amundson, D.E. 1993. Lack of agreement between measurement of ejection fraction by impedance cardiography versus radionuclide ventriculography. Critical Care Medicine 21: 1523–1527.


%% Copyright
% Copyright (C) 2017 Erik Huizinga, huizinga.erik@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Create fresh environment
%TODO deze sectie verwijderen
clc
clear
close all


%% Display demo text
disp 'Demonstration of Bland-Altman-Analysis'


%% Load data
% Add data folder to path. The following line only works when running this
% script or when the current directory is contains the ba19999data folder,
% otherwise a warning is issued.
pth = fileparts(mfilename('fullpath'));
addpath(fullfile(pth, 'ba1999data'))

% Load Table 1 [BA1999, p. 137-8]
load bpdata  % Systolic blood pressure (mmHg) data

% Extract data
n = BPData(:, 1);  % Subject number
J = BPData(:, 2 : 4);  % Three successive measurements by expert J
R = BPData(:, 5 : 7);  % Three successive measurements by expert R
S = BPData(:, 8 : 10);  % Three successive measurements by machine S
JName = 'Systolic BP, Observer J (mmHg)';
SName = 'Systolic BP, Machine S (mmHg)';


%% 2 Limits of agreement [BA1999, p. 136]
disp 'Section 2 Limits of agreement'

% Consider only the first observations by J and S
J1 = J(:, 1);
S1 = S(:, 1);

% Calculate Bland-Altman statistics
stats = ba(J1, S1);


% Show results, compare with comments
bias = stats.difference.bias;  % Bias (mean difference)
display(bias)  % [BA1999]: -16.29 mmHg

std = stats.difference.std;  % Standard deviation of difference
display(std)  % [BA1999]: 19.61 mmHg

loa = stats.difference.loa;  % Limits of agreement (LOA)
display(loa)  % [BA1999]: [-54.7, 22.1] mmHg


% Exclude subjects 78 and 80 [BA1999, p. 139]
% The indices of the outliers can be found using the data cursor: click on
% a sample in the graph to see the subject number i
fprintf '\nExcluding subjects 78 and 80 (outliers):\n'
iEx = [78, 80];
lEx = (n == 78) | (n == 80);  % Alternative: logical indices
stats = ba(J1, S1, 'Exclude', iEx);


% Show results, compare with comments
bias = stats.difference.bias;
display(bias)  % [BA1999]: -14.9 mmHg
% Erratum: a small error, probably a typo, exist in the article [BA1999].
% Note the difference between the bias and the article's mean difference.
% However, because the limits of agreement are the same and, of course, the
% mean difference equals the mean of the limits of agreement, thus the
% calculation here appears to be correct.

loa = stats.difference.loa;
display(loa)  % [BA1999]: [-43.6, 15.0] mmHg


%% 2.1 Graphical presentation of agreement [BA1999 p. 140]
disp 'Section 2.1 Graphical presentation of agreement'

% If any other figures exist, the numbers below might not be correct. This
% script was written with all other figures closed.
% Figure 1 corresponds to figure 1 in [BA1999], thus figure 2 corresponds
% to figure 2 in [BA1999].
f1_2 = figures(2);

% Perform Bland-Altman Analysis and create both the correlation and
% mean-difference graphs
stats = ba(f1_2, J1, S1, 'XName', JName, 'YName', SName, ...
           'PlotDefault', true);


% Show results, compare with comments
r = stats.difference.Spearman.r;
% [BA1999]: 0.07 (p. 140), here -0.03, so an error exists either here or in
% the article. Either way, the correlation is not significantly different
% from zero: see stats.difference.Spearman.p for the p-value.
display(r) % Spearman rank correlation between mean and difference


% Add limits of agreement
% This figure is equal to the previous, but has added limits of agreement
% It corresponds to figure 3 in [BA1999]
f3 = figure;

% Create the mean-difference graph and plot basic statistics, i.e. bias and
% limits of agreement
ba(f3, J1, S1, 'XName', JName, 'YName', SName, ...
    'PlotMeanDifference', true, 'PlotStatistics', 'basic')

% Print some text output
f = [f1_2; f3];
fn = sort([f.Number]);
disp(['See also figures [', num2str(fn), '].'])


%% 2.2 Precision of the estimated limits of agreement [BA1999, p. 141]
disp 'Section 2.2 Precision of the estimated limits of agreement'

% Clear figure 3 and recreate it with additional confidence intervals
% (not shown in article)
clf  

% Perform Bland-Altman analysis
stats = ba(f3, J1, S1, 'XName', JName, 'YName', SName, ...
           'PlotMeanDifference', true, 'PlotStatistics', 'extended');


% Show results, compare with comments
biasCI = stats.difference.biasCI;
% muDCI [BA1999, p. 142]: [-20.5 -12.1]
display(biasCI)

loaCI = stats.difference.loaCI; % confidence interval of the loa
% loaCI [BA1999, p. 142]:
% [-61.9, 14.9
%  -47.5, 29.3] mmHg
display(loaCI)

% Print some text output
f = f3;
fn = f.Number;
disp(['See also figure ', num2str(fn), '.'])


%% 3 Relationship between difference and magnitude [BA1999, p. 142]
disp 'Section 3 Relationship between difference and magnitude'

% Load Table 2 [BA1999, p. 143]
load pvdata  % Plasma volume in (%) data

% Extract data
% n = PVData(:, 1);
Nadler = PVData(:, 2);
Hurley = PVData(:, 3);
NName = 'Plasma volume (Nadler) (%)';
HName = 'Plasma volume (Hurley) (%)';

% Prepare figure and axes
% The following figure corresponds to figure 4 in [BA1999]. Its subplots
% correspond to figures 4a and 4b respectively.
f4 = figure;
ax4(1) = subplot(2, 1, 1);
ax4(2) = subplot(2, 1, 2);

% Perform Bland-Altman analysis
% The following could be accomplished with one line of code, but the graphs
% in the article are constructed a bit differently, so here we create the
% graphs accordingly, but need to do it with two separate lines.
ba(ax4(1), Hurley, Nadler, 'XName', HName, 'YName', NName, ...
    'PlotCorrelation', true)
ba(ax4(2), Nadler, Hurley, 'XName', NName, 'YName', HName, ...
    'PlotMeanDifference', true)
% In one line of code it would have been:
% ba(ax, Hurley, Nadler, 'PlotDefault', true)
% However, figure 4b in [BA1999] has the vertical axis flipped, so this
% does not work. Note that there is significant positive Spearman rank
% correlation between difference and mean, which can be seen in the legend
% if the 'PlotStatistics', 'basic' Name-Value pair argument is passed to
% ba, or by looking at the statistics in s, the output argument of s =
% ba(...).

% Print some text output
fn = f4.Number;
disp(['See also figure ', num2str(fn), '.'])


%% 3.1 Logarithmic transformation [BA1999, p. 143]
disp 'Section 3.1 Logarithmic transformation'

% Prepare figure and axes
% Figure 5 in [BA1999] is reproduced in the following figure. The subplots
% are article figures 5a and 5b respectively.
f5 = figure;
ax5(1) = subplot(2, 1, 1);
ax5(2) = subplot(2, 1, 2);

% Perform BAA using log transformation
ba(ax5(1), Hurley, Nadler, 'XName', HName, 'YName', NName, ...
    'PlotCorrelation', true, 'Transform', @log)
stats = ba(ax5(2), Nadler, Hurley, 'XName', NName, 'YName', HName, ...
           'PlotMeanDifference', true, ...
           'PlotStatistics', 'basic', ...
           'Transform', @log);


% Show results, compare with comments
logBias = stats.difference.bias;
display(logBias)  % [BA1999]: 0.099

logLoa = stats.difference.loa;
display(logLoa)  % [BA1999]: [0.056, 0.141]

logLoaCI = stats.difference.loaCI;
lowerLogLoaCI = logLoaCI(:, 1);  % Lower LOA CI
display(lowerLogLoaCI)  % [BA1999] [0.049; 0.064]

% Transform results back to regular units, compare with comments
bias = exp(logBias);
display(bias)  % [BA1999]: 1.11

loa = exp(logLoa);
display(loa)  % [BA1999]: [1.06, 1.15]


% Figure 6 in [BA1999] corresponds to this figure.
f6 = figure;

% Create the mean-ratio graph
ba(f6, Nadler, Hurley, 'XName', NName, 'YName', HName, ...
    'PlotMeanRatio', true, 'PlotStatistics', 'basic');

% Print some text output
f = [f5; f6];
fn = sort([f.Number]);
disp(['See also figures [', num2str(fn), '].'])


%% 3.2 A regression approach for nonuniform differences (article p. 145)
disp 'Section 3.2 A regression approach for nonuniform differences'

% Load Table 3 [BA1999, p. 146]
load fcdata

% Extract data
Trig = FCData(:, 1);
Gerber = FCData(:, 2);
TName = 'Fat (g/100 ml; Trig.)';
GName = 'Fat (g/100 ml; Gerber)';


% Prepare figure and axes
% Figure 7 [BA1999, p. 147] and its subfigures correspond to the following
% figure and its subplots.
f7 = figure;
ax7(1) = subplot(2, 1, 1);
ax7(2) = subplot(2, 1, 2);

% Perform default BAA on the data
ba(ax7(1), Gerber, Trig, 'XName', GName, 'YName', TName, ...
    'PlotCorrelation', true)
ba(ax7(2), Trig, Gerber, 'XName', TName, 'YName', GName, ...
    'PlotMeanDifference', true)


% Figure 8 [BA1999, p. 148] corresponds to the following figure.
f8 = figure;

% Perform BAA, but now plot the simple linear regression line of the
% difference on the mean instead of a constant bias. The corresponding 95%
% limits of agreement are plotted as well. In the article, the residuals of
% this regression line are not assumed to be significantly associated with
% the mean, so the assumption is made the residual variance is constant.
% This is specified using the 'ConstantResidualVariance', true Name-Value
% pair argument. The result of this assumption is the 95% LOA being
% parallel to the bias regression line.
stats = ba(f8, Trig, Gerber, 'XName', TName, 'YName', GName, ...
           'PlotMeanDifference', true, 'PlotStatistics', 'regression', ...
           'ConstantResidualVariance', true);


% Show results, compare with comments
% standard deviation of residual of simple linear regression polynomial of
% difference on mean:
stde = stats.difference.poly.stde;
display(stde)  % [BA1999]: 0.08033 (s_d on p. 148)
% Erratum: [BA1999] shows a slightly different value. This is less than
% 1.2% error, but the difference is there.

% Print some text output
f = [f7; f8];
fn = [f.Number];
disp(['See also figures [', num2str(fn), '].'])


%% 4 The importance of repeatability [BA1999, p. 148]
%TODO
% disp 'Section 4 The importance of repeatability'


%% 4.1 Estimating repeatability [BA1999, p. 149]
%TODO
% disp 'Section 4.1 Estimating repeatability'


%% 5 Measuring agreement using repeated measurements [BA1999, p. 150]
disp 'Section 5 Measuring agreement using repeated measurements'


%% 5.1 Equal numbers of replicates [BA1999, p. 150]
disp 'Section 5.1 Equal numbers of replicates'

% Perform Bland-Altman Analysis for repeated measurements. Notice the
% syntax is identical to the regular BAA, but because of the matrix input
% (check size(J) and size(S)) ba(J, S) performs BAA for repeated
% measurements.
stats = ba(J, S);


% Show results, compare with comments
varwJ = stats.xy.varw.x;  % Within-subject variance of measurements J
display(varwJ)  % [BA1999, p. 151]: 37.408

varwS = stats.xy.varw.y;  % Within-subject variance of measurements S
display(varwS)  % [BA1999, p. 151]: 83.141

bias = stats.difference.bias;  % Mean difference between J and S (bias)
display(bias)  % [BA1999, p. 151]: -15.62 mmHg

std = stats.difference.std;  % Standard deviation of the difference
display(std)  % [BA1999, p. 152]: 20.95 mmHg

loa = stats.difference.loa;  % Limits of agreement
display(loa)  % [BA1999, p. 152]: [-56.68, 25.44] mmHg
% Notice how the values of bias, std and loa are very similar to those of
% [BA1999] section 2, which is to be expected (they do not change for
% repeated measurements).

% Show more results, compare with comments
biasCI = stats.difference.biasCI;
display(biasCI)
% [BA1999]: value not presented, but note similarity to biasCI in section 2

loaCI = stats.difference.loaCI;
display(loaCI)  % [BA1999, p. 153]:
% [-63.5, 18.70
%  -49.9, 32.2] mmHg
% [BA1999] uses an incorrect calculation, because it uses 1.96. This value
% is the 97.5 percentile of the normal distribution, whereas Student's
% t-distribution should have been used. This would have resulted in a
% slightly larger value than 1.96, hence the article's estimation of the
% confidence intervals of the LOA is too narrow.


%% 5.2 Unequal numbers of replicates [BA1999, p. 154]
disp 'Section 5.2 Unequal numbers of replicates'

% Load Table 4 [BA1999, p. 154]
% This data was used by both [BA1999] and [BA2007] from [BSOCA1993].
% In [BA1999] the table (Table 4) differs from Table 1 in [BA2007]. The
% corresponding figure in [BA2007] (Figure 1) does not display the same
% data as the table! Furthermore, the data in the [BA1999] table seems to
% correspond to the figure in [BA2007]. This leads to the conclusion the 
% table in [BA2007] is incorrect and the table from [BA1999] is correct.
% Thus, codata.mat contains the [BA1999] cardiac output data.
% Unfortunately, the original data source [BSOCA1993] doesn't tabulate the
% data in the article. Bland and Altman probably got it directly from the
% authors.
load codata

% Extract and prepare data
n = COData(:, 1);
RV = COData(:, 2);
IC = COData(:, 3);

un = unique(n);
un = un(:).';
for i = numel(un) : -1 : 1
    j = un(i);
    cellRV{i} = RV(n == j);
    cellIC{i} = IC(n == j);
end

RVName = 'Radionuclide ventriculography';
ICName = 'Impedance cardiography';

% Figure 9 in [BA1999, p. 155] corresponds to the following figure
f9 = figure;
ax9(1) = subplot(2, 1, 1);
ax9(2) = subplot(2, 1, 2);

% Graph subject mean against within-subject standard deviation for each
% method separately
ba(ax9, cellRV, cellIC, 'XName', RVName, 'YName', ICName, ...
    'PlotMeanSD', 'separate', 'PlotStatistics', 'basic')


% Figure 10 in [BA1999, p. 156] corresponds to the following figure.
f10 = figure;

% Perform Bland-Altman Analysis for repeated measurements with unequal
% replicates
stats = ba(f10, cellRV, cellIC, 'XName', RVName, 'YName', ICName, ...
          'PlotMeanDifference', true);


% Show results, compare with comments
varwRV = stats.xy.varw.x;  % Within-subject variance component from RV
display(varwRV)  % [BA1999]: 0.1072

varwIC = stats.xy.varw.y;  % Within-subject variance component from IC
display(varwIC)  % [BA1999]: 0.1379

% Standard deviation of differences between single observations
std = stats.difference.std;
display(std)  % [BA1999]: 1.0517

bias = stats.difference.bias;  % Mean difference, i.e. bias
display(bias)  % [BA1999]: 0.7092

loa = stats.difference.loa;  % Limits of agreement
display(loa)  % [BA1999]: loa = [-1.3521, 2.7705]

% Show more results
biasCI = stats.difference.biasCI;
display(biasCI)  % [BA1999]:
% Value not presented, but given here for completeness

loaCI = stats.difference.loaCI;
display(loaCI)  % [BA1999]:
% Value not presented, but given here for completeness

% Some text output
f = [f9; f10];
fn = [f.Number];
disp(['See also figures [', num2str(fn), '].'])


%% 5.3 Replicated data in pairs [BA1999, p. 156]
disp 'Section 5.3 Replicated data in pairs'
% Note 1/2:
% This section is implemented in ba.m based on the [BA2007], because the
% section in [BA1999] contains an error. The [BA2007] article corrects this
% error and provides a calculation example with the same data as [BA1999].
% The numbers calculated below are the numbers from [BA2007, §3, p. 575
% onwards].

% Note 2/2:
% In [BA2007] Figure 3 the vertical axis depicts the difference (and this
% is plotted as well, hence the figure looks a lot like [BA1999] Figure
% 10). However, the caption mentions the standard deviation is plotted!
% This is incorrect. To plot the difference standard deviation versus the
% mean specify the 'PlotMeanSD', 'difference' Name-Value pair.

% Perform BAA
stats = ba(cellRV, cellIC, 'XName', RVName, 'YName', ICName, ...
           'ConstantTrueValue', false);


% Show results, compare with comments
std = stats.difference.std;
display(std)  % [BA2007]: 0.99062408
% (to see this many digits in MATLAB, use: format long)

bias = stats.difference.bias;
display(bias)  % [BA2007]: 0.6021667

loa = stats.difference.loa;
display(loa)  % [BA2007]: [-1.3394565, 2.5437899]


%% Demonstrate other functionality
% Create honeycomb plot (2D histogram with hexagonal bins) instead of
% regular scatter plot
f_11_12 = figures(2);
stats = ba( f_11_12, ...
    J1, S1, 'XName', JName, 'YName', SName, ...
    'PlotDefault', true, ...
    'PlotHoneycomb', true, ...
    'PlotStatistics', 'extended');

% Get the handles of the honeycomb patch objects
honeyHandle1 = stats.graphics.xy;
honeyHandle2 = stats.graphics.difference;


% Adjust honeycomb plot 1
honeyHandle1.EdgeColor = 'none';


% Adjust honeycomb plot 2
% Create color map from white to blue
map = lines;
map = map(1, :);
map = interp1([1, 0], [map; 1, 1, 1], linspace(0, 1, 256));

% Set color map
honeyHandle2.EdgeColor = map(end, :);
honeyHandle2.LineWidth = 2;
colormap(map)

% Add color bar
bar = colorbar;
ylabel(bar, 'counts', 'Rotation', -90, 'VerticalAlignment', 'bottom')