function [doPlotMD, axMD, doPlotMR, axMR, MSDType, doPlotMSD1, axMSD1, ...
          doPlotMSD2, axMSD2, doPlotC, axC] ...
          = validatePlotArgs(PlotDefault, PlotMeanDifference, ...
                             PlotMeanRatio, PlotMeanSD, PlotCorrelation, h)

                         
% Set default axes variables
axMD = [];
axMR = [];
axMSD1 = [];
axMSD2 = [];
axC = [];


% Set doPlotMD and doPlotC
doPlotDefault = logical(PlotDefault);
if doPlotDefault
    doPlotMD = true;  % Do mean-difference plot
    doPlotC = true;  % Do correlation plot
    
else
    doPlotMD = logical(PlotMeanDifference);
    doPlotC = logical(PlotCorrelation);
end


% Validate and parse mean-standard deviation graphs
diffstr = {'difference','single','joint'};
ratiostr = {'ratio'};
separatestr = {'separate','both','input'};
nonestr = {'none'};
doPlotMR = logical(PlotMeanRatio);
MSDType = parseMeanSD(PlotMeanSD);

switch MSDType
    case diffstr
        doPlotMSD1 = true;
        doPlotMSD2 = false;
        MSDType = 'difference';
    case ratiostr
        doPlotMSD1 = true;
        doPlotMSD2 = false;
        MSDType = 'ratio';
    case separatestr
        doPlotMSD1 = true;
        doPlotMSD2 = true;
        MSDType = 'separate';
    otherwise
        doPlotMSD1 = false;
        doPlotMSD2 = false;
end


% Validate number and type of handles in h for the requested plots
doPlot = [doPlotMD, doPlotMR, doPlotMSD1, doPlotMSD2, doPlotC];

if any(doPlot)
    nPlot = nnz(doPlot);
    if isempty(h)
        figure
        ax = gobjects(0);
        
    else
        nh = numel(h);
        if isaxes(h)
            if nh == nPlot
                ax = h;
                
            else
                error(['The number of requested plots (' num2str(nPlot) ...
                      ') is not equal to the number of axes (', ...
                      num2str(nh), ') in the first input argument.'])
            end
            
        else  % h must be a figure or contain nPlot figures
            if nh == 1
                figure(h)  % Throw error if h is not a figure
                ax = gobjects(0);
                
            elseif nh == nPlot
                for n = nh:-1:1
                    figure(h(n))
                    ax(n) = axes;
                end
                
            elseif all(isfigure(h))
                error(['The number of requested plots (' num2str(nPlot) ...
                      ') is not equal to the number of figures (', ...
                      num2str(nh), ') in the first input argument.'])
                
            else
                error(['The first optional input argument must be an ' ...
                      '(array of) axes or figure handle(s).'])
            end
        end
    end
    
    % Create axes if not yet existent
    if isempty(ax)
        if nPlot>1
            for n = nPlot:-1:1
                ax(n) = subplot(nPlot, 1, n);
            end
            
        else
            ax = axes;
        end
    end
    
    % store axes in plot order
    if doPlotMD, axMD = ax(1); ax(1) = []; end
    if doPlotMR, axMR = ax(1); ax(1) = []; end
    if doPlotMSD1, axMSD1 = ax(1); ax(1) = []; end
    if doPlotMSD2, axMSD2 = ax(1); ax(1) = []; end
    if doPlotC, axC = ax(1); end % ax(1) = []; end
end


    function str = parseMeanSD(in)
        % Parse PlotMeanSD
        if islogical(in)
            if in
                str = 'difference';
                
            else
                str = 'none';
            end
            
        else  % in must be a string
            % Validate
            str = validatestring(lower(in), ...
                [diffstr, ratiostr, separatestr, nonestr] ...
                );
        end
    end
end