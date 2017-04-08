%TESTSCATTERSET Unit tests for the ScatterSet Name-Value pair argument


%% Initialize
clear
close all


%% Load data
% Add data to path
pth = fileparts(mfilename('fullpath'));
firstChar = pth(1);
pth = strsplit(pth, filesep);
pth = pth(1 : end-1);
pth = fullfile(pth{:}, 'ba1999data');
if pth(1) ~= firstChar
    pth = [firstChar, pth];
end
addpath(pth)

% Load data
load bpdata

% Extract data
n = BPData(:, 1);  % Subject number
J = BPData(:, 2 : 4);  % Three successive measurements by expert J
R = BPData(:, 5 : 7);  % Three successive measurements by expert R
S = BPData(:, 8 : 10);  % Three successive measurements by machine S
JName = 'Systolic BP, Observer J (mmHg)';
SName = 'Systolic BP, Machine S (mmHg)';


%% Test plotting data without repeated measurements
% Consider only the first observations by J and S
J1 = J(:, 1);
S1 = S(:, 1);

% Test full set: nothing should happen
f = figure;
ba(f, J1, S1, 'XName', 'J1, no reps full', 'YName', 'S1, no reps full', ...
    'PlotDefault', true, 'PlotStatistics', 'extended', ...
    'ScatterSet', 'full')
fprintf(['Figure ', num2str(f.Number), ' shows regular BAA with ', ...
         'the full scatter set,\n  i.e., all observations, i.e., ', ...
         'one observation per subject.\n'])


% Test mean set: nothing should happen, because the mean of a measurement
% without repetitions equals the measurement itself
f = figure;
ba(f, J1, S1, 'XName', 'J1, no reps mean', 'YName', 'S1, no reps mean', ...
    'PlotDefault', true, 'PlotStatistics', 'extended', ...
    'ScatterSet', 'mean')
fprintf(['Figure ', num2str(f.Number), ' shows regular BAA with ', ...
         'the mean scatter set,\n  i.e., all observations, i.e., ', ...
         'one observation per subject.\n'])


%% Test plotting full set
f = figure;
ba(f, J, S, 'XName', 'J, full', 'YName', 'S, full', ...
    'PlotDefault', true, 'PlotStatistics', 'extended', ...
    'ScatterSet', 'full')
fprintf(['Figure ', num2str(f.Number), ' shows BAA for repeated ', ...
         'measurements with the full scatter set,\n  i.e., all ', ...
         'observations, i.e., multiple observations per subject.\n'])


%% Test plotting means in set
f = figure;
ba(f, J, S, 'XName', 'J, mean', 'YName', 'S, mean', ...
    'PlotDefault', true, 'PlotStatistics', 'extended', ...
    'ScatterSet', 'mean')
fprintf(['Figure ', num2str(f.Number), ' shows BAA for repeated ', ...
         'measurements with the mean scatter set,\n  i.e., one mean ', ...
         'per subject.\n'])