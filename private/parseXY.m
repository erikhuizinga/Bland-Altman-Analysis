function [xok, yok, doRepeated] = parseXY(x, y)

% Check validity of x and y
if iscell(x) && iscell(y)
    lok = cellfun(@(x,y) isfinite(x) & isnumeric(x) ...
                  & isfinite(y) & isnumeric(y), ...
                  x, y, 'UniformOutput',false);
    
else
    lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
end


% Parse x and y into arrays of valid observations, determine repeated
% measurements
[xok, repX] = parseV(x, lok);
[yok, repY] = parseV(y, lok);
doRepeated = repX | repY;


% Check number of elements
if isvector(xok) && isvector(yok)
    % x and y are numeric or cells
    assert(numel(xok) == numel(yok), ...
           'Number of elements in x and y must be equal.')
       
else  % x and/or y are matrices
    % Check number of rows (subjects)
    assert(size(xok, 1) == size(yok, 1), ...
           ['x and y must have the same number of rows (subjects) for ' ...
           repeated measures.'])
    
    % Check number of columns (replicates)
    assert(size(xok, 2) == size(yok, 2), ...
           ['x and y must have the same number of columns ' ...
           '(replicates) per subject.'])
end
end


function [vok, doRepeated] = parseV(v, lok)
if iscell(v)
    doRepeated = true;
    
    
    % Check cell contents
    assert(all(cellfun(@isvector, v)), ...
           ['For cell input, every cell of ', inputname(1), ...
           ' must contain a vector.'])
    
    
    % Reshape cell contents
    vok = cellfun(@(v) v(:).', v, 'UniformOutput', false);
    vok = vok(:);
    % Now vok is a column vector cell, every element a subject and every
    % cell's contents a row vector of replicates
    
    
    % Keep valid indices
    % The following apparently doesn't ‘reshape’:
    % reshape(lok,size(vok));  % Match shapes
    if ~isequal(size(lok), size(vok)), lok = lok.'; end  % Match shapes
    vok = cellfun(@(vok, lok) vok(lok), vok, lok, 'UniformOutput', false);
    
    
    % Determine number of valid elements in cells
    nV = cellfun(@numel, vok);
    
    if isscalar(unique(nV))
        % Reshape to a matrix, because all cells contain the same number of
        % replicates
        vok = vertcat(vok{:});
        % Now v is a numeric matrix with rows corresponding to subjects and
        % columns corresponding to replicates
        
        if isvector(vok)
            % Reparse, because vok contains a single measurement per cell
            % element, thus no replicates
            [vok, doRepeated] = parseV(vok, true(size(vok)));
        end
    end
    
else  % v is not a cell
    if isvector(v)
        vok = v(lok);
        vok = vok(:);  % Force column vectors
        doRepeated = false;
        
    else  % v is a numeric matrix
        if isscalar(unique(sum(lok, 2)))
            vok = v(lok);
            nCol = size(lok, 2);
            vok = reshape(vok, [], nCol);
            % It is assumed the user provides repeated observations, every
            % row being a subject and every column an observation
            
        else
            % Reshape x and y into cells, every row a subject and every
            % cell a vector of replicates
            for r = flip(find(any(lok, 2))).'
                vok{r, 1} = v(r, lok(r, :));  %#ok<AGROW>
            end
        end
        
        doRepeated = true;
    end
end
end