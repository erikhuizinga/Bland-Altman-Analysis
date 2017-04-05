function lineColor = getLineColor(name)
switch lower(name)
    case 'scatter'
        lineColor = [0, 0, 0];  % Black
        
    case 'honeycomb'
        lineColor = [1, 0, 0];  % Red
        
    otherwise
        error 'Unknown scatter function name'
end