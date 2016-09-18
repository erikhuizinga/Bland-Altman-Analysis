%% load data
% add data folder to path
pth = fileparts(mfilename('fullpath'));
addpath(fullfile(pth,'ba1999data'))

% load Table 1 (article p. 137-8)
load bpdata % systolic blood pressure (mmHG) data

n = BPData(:,1); % subject number
J = BPData(:,2:4); % three successive measurements by expert J
R = BPData(:,5:7); % three successive measurements by expert R
S = BPData(:,8:10); % three successive measurements by machine S
JName = 'Systolic BP, Observer J (mmHg)';
SName = 'Systolic BP, Machine S (mmHg)';

%% 5 Measuring agreement using repeated measurements (article p. 150)
disp 'Section 5 Measuring agreement using repeated measurements'
%% 5.1 Equal numbers of replicates (article p. 150)
disp 'Section 5.1 Equal numbers of replicates'

% perform Bland-Altman Analysis for repeated measurements. Notice the
% syntax is identical to the regular BAA, but because of the matrix input
% (check size(J) and size(S)) ba(J,S) performs BAA for repeated
% measurements.
s = ba(J,S);

% show results
varWithinJ = s.x.varWithin; % within-subject variance of measurements J
varWithinS = s.y.varWithin; % within-subject variance of measurements S
display(varWithinJ) % BA1999 p. 151: 37.408
display(varWithinS) % BA1999 p. 151: 83.141
muD = s.difference.mu; % mean difference between J and S (bias)
display(muD) % BA1999 p. 151: -15.62 mmHg
sD = s.difference.s; % standard deviation of the difference
display(sD) % BA1999 p. 152: 20.95 mmHg
loaD = s.difference.loa; % limits of agreement
display(loaD); % BA1999 p. 152: [-56.68, 25.44] mmHg
% Notice how the values of muD, sD and loa are very similar to those of
% section 2, which is to be expected (they do not change for repeated
% measurements).

% more results
muDCI = s.difference.muCI;
display(muDCI); % BA1999:
% value not presented, but notice similarity with muDCI in section 2.
loaDCI = s.difference.loaCI;
display(loaDCI); % BA1999 p. 153:
% [-63.5, 18.70
%  -49.9, 32.2] mmHg
% The article uses an incorrect calculation, because it uses 1.96. This
% value comes from the normal distribution, whereas Student's
% t-distribution should have been used. This would have resulted in a
% slightly larger value than 1.96, hence the article's estimation of the
% confidence intervals of the LOA is too narrow.

%% 5.2 Unequal numbers of replicates (article p. 154)
disp 'Section 5.2 Unequal numbers of replicates'

% load Table 4 (article p. 154)
load codata
n = COData(:,1);
RV = COData(:,2);
IC = FCData(:,3);
RVName = 'Radionuclide ventriculography';
ICName = 'Impedance cardiography';

%% 5.3 Replicated data in pairs (article p. 156)
disp 'Section 5.3 Replicated data in pairs'