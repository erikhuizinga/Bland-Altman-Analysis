function [doPlotBasicStats, doPlotExtendedStats, doPlotRegStats, ...
          doConstantRegression] ...
          = parseStatArgs(PlotStatistics, ConstantResidualVariance)


% Parse statistics arguments

% Set defaults
doPlotBasicStats = false;
doPlotExtendedStats = false;
doPlotRegStats = false;
doConstantRegression = logical(ConstantResidualVariance);

switch lower(PlotStatistics)
    case 'none'
        % Keep defaults
    case 'basic'
        doPlotBasicStats = true;
    case 'extended'
        doPlotBasicStats = true;
        doPlotExtendedStats = true;
    case 'regression'
        doPlotBasicStats = true;
        doPlotRegStats = true;
    otherwise
        error 'Unknown value for the ''PlotStatistics'' name-value pair.'
end