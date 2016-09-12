%BA1999DEMO Demonstration of:
% Bland and Altman 1999
% Measuring agreement methods in comparison studies
% http://smm.sagepub.com/content/8/2/135.abstract

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

%% display demo text
disp 'Demonstration of Bland-Altman analysis'

%% load data
% add data folder to path
pth = fileparts(mfilename('fullpath'));
addpath(fullfile(pth,'ba1999data'))

% load Table 1 (article p. 137-8)
load bpdata

n = BPData(:,1); % subject number
J = BPData(:,2:4); % three successive measurements by expert J
R = BPData(:,5:7); % three successive measurements by expert R
S = BPData(:,8:10); % three successive measurements by machine S
JName = 'Systolic BP, Observer J (mmHg)';
SName = 'Systolic BP, Machine S (mmHg)';

%% 2 Limits of agreement (article page 136)
disp 'Section 2 Limits of agreement'

% This section only considers the first observations by J and S.
J1 = J(:,1);
S1 = S(:,1);

s = ba(J1,S1);
muD = s.muD; % mean difference (bias)
sD = s.sD; % standard deviation of difference
loa = s.loa; % limits of agreement
display(muD) % article: -16.29 mmHg
display(sD) % article: 19.61 mmHg
display(loa) % article: [-54.7, 22.1] mmHg

% exclude subjects 78 and 80 (article page number 139)
fprintf '\nExcluding subjects 78 and 80 (outliers):\n'
iEx = [78,80];
lEx = (n==78) | (n==80); % alternative: logical indices
s = ba(J1,S1, 'Exclude',iEx);
muD = s.muD;
loa = s.loa;
display(muD) % article: -14.9 mmHg
% Erratum: a small error, probably a typo, exist in the original article.
% Note the difference between the muD and the article's mean difference.
% However, because the limits of agreement are the same and, of course, the
% mean difference equals the mean of the limits of agreement, thus the
% calculation here appears to be correct.
display(loa) % article: [-43.6, 15.0] mmHg

%% 2.1 Graphical presentation of agreement (article page 140)
disp 'Section 2.1 Graphical presentation of agreement'

% If any other figures exist, the numbers below might not be correct. This
% script was written with all other figures closed.
% Figure 1 corresponds to article figure 1, figure 2 corresponds to article
% figure 2.
f = figures(2);
s = ba(f, J1,S1, 'XName',JName, 'YName',SName, 'PlotAll',true);
rSMuD = s.rSMuD;
% Article: 0.07 (p. 140), here -0.03, so an error exists either here or in
% the article.
display(rSMuD); % Spearman rank correlation between mean and difference

% Figure 3 is equal to figure 2, but has added limits of agreement. It
% corresponds to article figure 3.
ba(figure, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','basic')

% some text output
f = [f;gcf];
fn = sort([f.Number]);
disp(['See also figures [' num2str(fn) '].'])

%% 2.2 Precision of the estimated limits of agreement (article page 141)
disp 'Section 2.2 Precision of the estimated limits of agreement'

clf % Clear figure 3 and recreate it with additional confidence intervals
% (not shown in article).
s = ba(gcf, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','extended');
muDCI = s.muDCI;
loaCI = s.loaCI; % confidence interval of the loa
% article loaCI (p. 142):
% [-61.9, 14.9
%  -47.5, 29.3]
display(loaCI)
% article muDCI (p. 142): [-20.5 -12.1]
display(muDCI)

% some text output
f = gcf;
fn = f.Number;
disp(['See also figure ' num2str(fn) '.'])

%% 3 Relationship between difference and magnitude (article p. 142)
disp 'Section 3 Relationship between difference and magnitude'

% load data from Table 2 in the article
load pvdata

n = PVData(:,1);
Nadler = PVData(:,2);
Hurley = PVData(:,3);

% create figures 4a and 4b in the article
f = figure;
ax(1) = subplot(1,2,1);
ax(2) = subplot(1,2,2);
% The following could be accomplished with one line of code, but the
% graphs in the article are constructed a bit differently, so here we
% create the graphs accordingly, but need to do it with two separate lines.
ba(ax(1), Hurley,Nadler, ...
    'XName','Plasma volume (Hurley) (%)', ...
    'YName','Plasma volume (Nadler) (%)', ...
    'PlotCorrelation',true)
ba(ax(2), Nadler,Hurley, ...
    'XName','Plasma volume (Nadler) (%)', ...
    'YName','Plasma volume (Hurley) (%)', ...
    'PlotMeanDifference',true)
% In one line of code it would have been:
% ba(ax, Hurley,Nadler, 'PlotAll',true)
% But figure 4b has the vertical axis flipped, so this does not work.

% some text output
fn = f.Number;
disp(['See also figure ' num2str(fn) '.'])

%% 3.1 Logarithmic transformation (article p. 143)
disp 'Section 3.1 Logarithmic transformation'

% perform BAA using log transformation
f = figure;
ax(1) = subplot(1,2,1);
ax(2) = subplot(1,2,2);
ba(ax(1), Hurley,Nadler, ...
    'XName','Log plasma volume (Hurley)', ...
    'YName','Log plasma volume (Nadler)', ...
    'PlotCorrelation',true, 'Transform',@log)
s = ba(ax(2), Nadler,Hurley, ...
    'XName','Log plasma volume (Hurley)', ...
    'YName','Log plasma volume (Nadler)', ...
    'PlotMeanDifference',true, ...
    'PlotStatistics','basic', ...
    'Transform',@log);
logMuD = s.muD;
display(logMuD) % article: 0.099
logLoa = s.loa;
display(logLoa) % article: [0.056, 0.141]
logLoaCI = s.loaCI;
lowerLogLoaCI = logLoaCI(:,1); % lower LOA CI
display(lowerLogLoaCI) % article [0.049, 0.064]

% backtransformation of results
muD = exp(logMuD);
display(muD) % article: 1.11
loa = exp(logLoa);
display(loa) % article: [1.06, 1.15]

% ratio mean-difference plot
ba(figure, Nadler,Hurley, ...
    'XName','Plasma volume (Nadler) (%)', ...
    'YName','Plasma volume (Hurley) (%)', ...
    'PlotMeanDifference',true, 'PlotStatistics','basic', ...
    'Transform','ratio' );

% some text output
f = [f;gcf];
fn = sort([f.Number]);
disp(['See also figures [' num2str(fn) '].'])

%% 3.2 A regression approach for nonuniform differences (article p. 145)
disp 'Section 3.2 A regression approach for nonuniform differences'