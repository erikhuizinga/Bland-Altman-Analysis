%BADEMO Demonstration of Bland-Altman-Analysis, a MATLAB implementation by
%   Erik Huizinga based on the methods by Bland and Altman 1983, 1986, 1999
%   and 2007. The references below are used throughout this demo script.
%
%% References
% [BA1999] Bland, J.M., and Altman, D.G. 1999. Measuring agreement in method comparison studies. Statistical Methods in Medical Research 8: 135-160. doi:10.1191/096228099673819272.
% [BA2007] Bland, J.M., and Altman, D.G. 2007. Agreement Between Methods of Measurement with Multiple Observations Per Individual. Journal of Biopharmaceutical Statistics 17: 571-582. doi:10.1080/10543400701329422.
% [BSOCA1993] Bowling, L.S., Sageman, W.S., O’Connor, S.M., Cole, R., and Amundson, D.E. 1993. Lack of agreement between measurement of ejection fraction by impedance cardiography versus radionuclide ventriculography. Critical Care Medicine 21: 1523–1527.

%% Copyright
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

%% fresh environment
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
addpath(fullfile(pth,'ba1999data'))

% load Table 1 [BA1999, p. 137-8]
load bpdata % systolic blood pressure (mmHG) data

n = BPData(:,1); % subject number
J = BPData(:,2:4); % three successive measurements by expert J
R = BPData(:,5:7); % three successive measurements by expert R
S = BPData(:,8:10); % three successive measurements by machine S
JName = 'Systolic BP, Observer J (mmHg)';
SName = 'Systolic BP, Machine S (mmHg)';

%% 2 Limits of agreement [BA1999, p. 136]
disp 'Section 2 Limits of agreement'

% This section only considers the first observations by J and S.
J1 = J(:,1);
S1 = S(:,1);

% calculate Bland-Altman statistics
stats = ba(J1,S1);

% show results, compare with comments
muD = stats.difference.mu; % mean difference (bias)
display(muD) % [BA1999]: -16.29 mmHg

sD = stats.difference.s; % standard deviation of difference
display(sD) % [BA1999]: 19.61 mmHg

loaD = stats.difference.loa; % limits of agreement (LOA)
display(loaD) % [BA1999]: [-54.7, 22.1] mmHg

% exclude subjects 78 and 80 [BA1999, p. 139]
% The indices of the outliers can be found using the data cursor: click on
% a sample in the graph to see the subject number i.
fprintf '\nExcluding subjects 78 and 80 (outliers):\n'
iEx = [78,80];
lEx = (n==78) | (n==80); % alternative: logical indices
stats = ba(J1,S1, 'Exclude',iEx);

% show results, compare with comments
muD = stats.difference.mu;
display(muD) % [BA1999]: -14.9 mmHg
% Erratum: a small error, probably a typo, exist in the article [BA1999].
% Note the difference between the muD and the article's mean difference.
% However, because the limits of agreement are the same and, of course, the
% mean difference equals the mean of the limits of agreement, thus the
% calculation here appears to be correct.

loaD = stats.difference.loa;
display(loaD) % [BA1999]: [-43.6, 15.0] mmHg

%% 2.1 Graphical presentation of agreement [BA1999 p. 140]
disp 'Section 2.1 Graphical presentation of agreement'

% If any other figures exist, the numbers below might not be correct. This
% script was written with all other figures closed.
% Figure 1 corresponds to figure 1 in [BA1999], thus figure 2 corresponds
% to figure 2 in [BA1999].
f1_2 = figures(2);

% Perform Bland-Altman Analysis and create both the correlation and
% mean-difference graphs.
stats = ba(f1_2, J1,S1, 'XName',JName, 'YName',SName, 'PlotDefault',true);

% show results, compare with comments
rSMuD = stats.difference.rSMu;
% [BA1999]: 0.07 (p. 140), here -0.03, so an error exists either here or in
% the article. Either way, the correlation is not significantly different
% from zero: see s.difference.pRSMu for the p-value.
display(rSMuD) % Spearman rank correlation between mean and difference

% This figure is equal to the previous, but has added limits of agreement.
% It corresponds to figure 3 in [BA1999].
f3 = figure;

% Create the mean-difference graph and plot basic statistics, i.e. bias and
% limits of agreement.
ba(f3, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','basic')

% some text output
f = [f1_2;f3];
fn = sort([f.Number]);
disp(['See also figures [' num2str(fn) '].'])

%% 2.2 Precision of the estimated limits of agreement [BA1999, p. 141]
disp 'Section 2.2 Precision of the estimated limits of agreement'

clf % Clear figure 3 and recreate it with additional confidence intervals
% (not shown in article).
stats = ba(f3, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','extended');

% show results, compare with comments
muDCI = stats.difference.muCI;
% muDCI [BA1999, p. 142]: [-20.5 -12.1]
display(muDCI)

loaDCI = stats.difference.loaCI; % confidence interval of the loa
% loaCI [BA1999, p. 142]:
% [-61.9, 14.9
%  -47.5, 29.3] mmHg
display(loaDCI)

% some text output
f = f3;
fn = f.Number;
disp(['See also figure ' num2str(fn) '.'])

%% 3 Relationship between difference and magnitude [BA1999, p. 142]
disp 'Section 3 Relationship between difference and magnitude'

% load Table 2 [BA1999, p. 143]
load pvdata % plasma volume in (%) data

% n = PVData(:,1);
Nadler = PVData(:,2);
Hurley = PVData(:,3);
NName = 'Plasma volume (Nadler) (%)';
HName = 'Plasma volume (Hurley) (%)';

% The following figure corresponds to figure 4 in [BA1999]. Its subplots
% correspond to figures 4a and 4b respectively.
f4 = figure;
ax(1) = subplot(2,1,1);
ax(2) = subplot(2,1,2);

% The following could be accomplished with one line of code, but the
% graphs in the article are constructed a bit differently, so here we
% create the graphs accordingly, but need to do it with two separate lines.
ba(ax(1), Hurley,Nadler, 'XName',HName, 'YName',NName, ...
    'PlotCorrelation',true)
ba(ax(2), Nadler,Hurley, 'XName',NName, 'YName',HName, ...
    'PlotMeanDifference',true)
% In one line of code it would have been:
% ba(ax, Hurley,Nadler, 'PlotAll',true)
% But figure 4b in [BA1999] has the vertical axis flipped, so this does not
% work. Note that there is significant positive Spearman rank correlation
% between difference and mean, which can be seen in the legend if the
% 'PlotStatistics','basic' Name-Value pair argument is passed to ba, or by
% looking at the statistics in s, the output argument of s = ba(...).

% some text output
fn = f4.Number;
disp(['See also figure ' num2str(fn) '.'])

%% 3.1 Logarithmic transformation [BA1999, p. 143]
disp 'Section 3.1 Logarithmic transformation'

% Figure 5 in [BA1999] is reproduced in the following figure. The subplots
% are article figures 5a and 5b respectively.
f5 = figure;
ax(1) = subplot(2,1,1);
ax(2) = subplot(2,1,2);

% perform BAA using log transformation
ba(ax(1), Hurley,Nadler, 'XName',HName, 'YName',NName, ...
    'PlotCorrelation',true, 'Transform',@log)
stats = ba(ax(2), Nadler,Hurley, 'XName',NName, 'YName',HName, ...
    'PlotMeanDifference',true, ...
    'PlotStatistics','basic', ...
    'Transform',@log);

% show results, compare with comments
logMuD = stats.difference.mu;
display(logMuD) % [BA1999]: 0.099

logLoaD = stats.difference.loa;
display(logLoaD) % [BA1999]: [0.056, 0.141]

logLoaDCI = stats.difference.loaCI;
lowerLogLoaDCI = logLoaDCI(:,1); % lower LOA CI
display(lowerLogLoaDCI) % [BA1999] [0.049; 0.064]

% backtransformation of results, compare with comments
muD = exp(logMuD);
display(muD) % [BA1999]: 1.11

loaD = exp(logLoaD);
display(loaD) % [BA1999]: [1.06, 1.15]

% Figure 6 in [BA1999] corresponds to this figure.
f6 = figure;

% Create the mean-ratio graph
ba(f6, Nadler,Hurley, 'XName',NName, 'YName',HName, ...
    'PlotMeanRatio',true, 'PlotStatistics','basic');

% some text output
f = [f5;f6];
fn = sort([f.Number]);
disp(['See also figures [' num2str(fn) '].'])

%% 3.2 A regression approach for nonuniform differences (article p. 145)
disp 'Section 3.2 A regression approach for nonuniform differences'

% load Table 3 [BA1999, p. 146]
load fcdata
Trig = FCData(:,1);
Gerber = FCData(:,2);
TName = 'Fat (g/100 ml; Trig.)';
GName = 'Fat (g/100 ml; Gerber)';

% Figure 7 [BA1999, p. 147] and its subfigures correspond to the following
% figure and its subplots.
f7 = figure;
ax(1) = subplot(2,1,1);
ax(2) = subplot(2,1,2);

% Perform default BAA on the data.
ba(ax(1), ...
    Gerber,Trig, 'XName',GName, 'YName',TName, ...
    'PlotCorrelation',true)
ba(ax(2), ...
    Trig,Gerber, 'XName',TName, 'YName',GName, ...
    'PlotMeanDifference',true)

% Figure 8 [BA1999, p. 148] corresponds to the following figure.
f8 = figure;

% Perform BAA, but now plot the simple linear regression line of the
% difference on the mean instead of a constant bias. The corresponding 95%
% limits of agreement are plotted as well. In the article, the residuals of
% this regression line are not assumed to be significantly associated with
% the mean, so the assumption is made the residual variance is constant.
% This is specified using the 'ConstantResidualVariance',true Name-Value
% pair argument. The result of this assumption is the 95% LOA being
% parallel to the bias regression line.
stats = ba(f8, Trig,Gerber, 'XName',TName, 'YName',GName, ...
    'PlotMeanDifference',true, ...
    'PlotStatistics','regression', 'ConstantResidualVariance',true);

% show results, compare with comments
% standard deviation of residual of simple linear regression polynomial of
% difference on mean:
sPolyResidualD = stats.difference.sPolyResidual;
display(sPolyResidualD) % [BA1999]: 0.08033 (s_d on p. 148)
% Erratum: [BA1999] shows a slightly different value. This is less than
% 1.2% error, but the difference is there.

% some text output
f = [f7;f8];
fn = [f.Number];
disp(['See also figures [' num2str(fn) '].'])

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
% (check size(J) and size(S)) ba(J,S) performs BAA for repeated
% measurements.
stats = ba(J,S);

% show results, compare with comments
varWithinJ = stats.x.varWithin; % within-subject variance of measurements J
display(varWithinJ) % [BA1999, p. 151]: 37.408

varWithinS = stats.y.varWithin; % within-subject variance of measurements S
display(varWithinS) % [BA1999, p. 151]: 83.141

muD = stats.difference.mu; % mean difference between J and S (bias)
display(muD) % [BA1999, p. 151]: -15.62 mmHg

sD = stats.difference.s; % standard deviation of the difference
display(sD) % [BA1999, p. 152]: 20.95 mmHg

loaD = stats.difference.loa; % limits of agreement
display(loaD) % [BA1999, p. 152]: [-56.68, 25.44] mmHg
% Notice how the values of muD, sD and loa are very similar to those of
% [BA1999] section 2, which is to be expected (they do not change for
% repeated measurements).

% show more results, compare with comments
muDCI = stats.difference.muCI;
display(muDCI) % [BA1999]:
% value not presented, but notice similarity to muDCI in section 2.

loaDCI = stats.difference.loaCI;
display(loaDCI) % [BA1999, p. 153]:
% [-63.5, 18.70
%  -49.9, 32.2] mmHg
% [BA1999] uses an incorrect calculation, because it uses 1.96. This value
% is the 97.5 percentile of the normal distribution, whereas Student's
% t-distribution should have been used. This would have resulted in a
% slightly larger value than 1.96, hence the article's estimation of the
% confidence intervals of the LOA is too narrow.

%% 5.2 Unequal numbers of replicates [BA1999, p. 154]
disp 'Section 5.2 Unequal numbers of replicates'

% load Table 4 [BA1999, p. 154]
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
n = COData(:,1);
RV = COData(:,2);
IC = COData(:,3);
un = unique(n);
un = un(:).';
for i = numel(un):-1:1
    j = un(i);
    cellRV{i} = RV(n==j);
    cellIC{i} = IC(n==j);
end
RVName = 'Radionuclide ventriculography';
ICName = 'Impedance cardiography';

% Figure 9 in [BA1999, p. 155] corresponds to the following figure.
f9 = figure;
ax(1) = subplot(2,1,1);
ax(2) = subplot(2,1,2);

% Graph subject mean against within-subject standard deviation for each
% method separately.
ba(ax, cellRV,cellIC, 'XName',RVName, 'YName',ICName, ...
    'PlotMeanSD','separate', 'PlotStatistics','basic')

% Figure 10 in [BA1999, p. 156] corresponds to the following figure.
f10 = figure;

% Perform Bland-Altman Analysis for repeated measurements with unequal
% replicates.
stats = ba(f10, cellRV,cellIC, 'XName',RVName, 'YName',ICName, ...
    'PlotMeanDifference',true);

% show results, compare with comments
varWithinRV = stats.x.varWithin; % within-subject variance component from RV
display(varWithinRV) % [BA1999]: 0.1072

varWithinIC = stats.y.varWithin; % within-subject variance component from IC
display(varWithinIC) % [BA1999]: 0.1379

% standard deviation of differences between single observations
sD = stats.difference.s;
display(sD) % [BA1999]: 1.0517

muD = stats.difference.mu; % mean difference, i.e. bias
display(muD) % [BA1999]: 0.7092

loa = stats.difference.loa; % limits of agreement
display(loa) % [BA1999]: loa = [-1.3521, 2.7705]

% show more results
muDCI = stats.difference.muCI;
display(muDCI) % [BA1999]:
% value not presented, but given here for completeness

loaDCI = stats.difference.loaCI;
display(loaDCI) % [BA1999]:
% value not presented, but given here for completeness

% some text output
f = [f9;f10];
fn = [f.Number];
disp(['See also figures [' num2str(fn) '].'])

%% 5.3 Replicated data in pairs [BA1999, p. 156]
disp 'Section 5.3 Replicated data in pairs'
% Note 1/2:
% This section is implemented in ba.m based on the [BA2007], because the
% section in [BA1999] contains an error. [BA2007] article corrects this
% error and provides a calculation example with the same data as [BA1999].
% The numbers calculated below are the numbers from [BA2007, §3, p. 575
% onwards].

% Note 2/2:
% In [BA2007] Figure 3 the vertical axis depicts the difference (and this
% is plotted as well, hence the figure looks a lot like [BA1999] Figure
% 10). However, the caption mentions the standard deviatino is plotted!
% This is incorrect. To plot the difference standard deviation versus the
% mean specify the 'PlotMeanSD','difference' Name-Value pair.

% perform BAA
stats = ba( cellRV,cellIC, 'XName',RVName, 'YName',ICName, ...
    'ConstantTrueValue',false);

% show results, compare with comments
sD = stats.difference.s;
display(sD) % [BA2007]: 0.99062408
% (to see this many digits in MATLAB, use format long)

muD = stats.difference.mu;
display(muD) % [BA2007]: 0.6021667

loaD = stats.difference.loa;
display(loaD) % [BA2007]: [-1.3394565, 2.5437899]