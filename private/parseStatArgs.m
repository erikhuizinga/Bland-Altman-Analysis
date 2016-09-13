function [doPlotBasicStats,doPlotExtStats,doPlotRegStats] = ...
    parseStatArgs(PlotStatistics)
% parse statistics arguments

% defaults
doPlotBasicStats = false;
doPlotExtStats = false;
doPlotRegStats = false;

switch lower(PlotStatistics)
    case 'none'
        % keep defaults
    case 'basic'
        doPlotBasicStats = true;
    case 'extended'
        doPlotBasicStats = true;
        doPlotExtStats = true;
    case 'regression'
        doPlotBasicStats = true;
        doPlotRegStats = true;
    otherwise
        error 'Unknown value for the ''PlotStatistics'' name-value pair.'
end