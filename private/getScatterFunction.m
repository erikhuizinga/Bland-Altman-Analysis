function scatterFunction = getScatterFunction(name)
import submodules.honeycomb.*

switch lower(name)
    case 'scatter'
        scatterFunction = @scatter;
        
    case 'honeycomb'
        scatterFunction = @honeycomb;
        
    otherwise
        error 'Unknown scatter function name'
end