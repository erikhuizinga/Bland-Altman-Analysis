function [xok,yok,doRepeated] = parseXY(x,y)

% check validity of x and y for calculations
if iscell(x) && iscell(y)
    lok = cellfun( ...
        @(x,y) isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y), ...
        x,y, 'UniformOutput',false);
else
    lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
end

[xok,repX] = parseV(x,lok);
[yok,repY] = parseV(y,lok);

doRepeated = repX|repY;

% check number of elements
if isvector(xok) && isvector(yok)
    % x and y are numeric or cells
    assert(numel(xok)==numel(yok), ...
        'Number of elements in x and y must be equal.')
else % x and/or y matrices
    % validate number of rows (subjects)
    assert(size(xok,1)==size(yok,1),['x and y must have the same ' ...
        'number of rows (subjects) for repeated measures.'])
    % check number of columns (replicates)
    assert(size(xok,2)==size(yok,2),['x and y must have the same ' ...
        'number of columns (replicates) per subject.'])
end
end

function [vok,doRepeated] = parseV(v,lok)
if iscell(v)
    doRepeated = true;
    
    assert(all(cellfun(@isvector,v)),['For cell input, every cell of ' ...
        inputname(1) ' must contain a vector.'])
    
    % reshape
    vok = cellfun(@(v) v(:).',v,'UniformOutput',false); % force row vectors
    vok = vok(:); % force column cell
    % now v is a column vector cell, every element a subject and every
    % cell's contents a row vector of replicates
    
    vok = cellfun(@(vok,lok) vok(lok), vok,lok, 'UniformOutput',false);
    
    nV = cellfun(@numel,vok);
    if isscalar(unique(nV))
        % all cells contain the same number of replicates, thus they fit in
        % a matrix
        vok = vertcat(vok{:}); % convert to matrix
        % now v is a numeric matrix with rows corresponding to subjects and
        % columns corresponding to replicates
        if isvector(vok)
            % a single measurement per cell element, thus no replicates,
            % reparse
            [vok,doRepeated] = parseV(vok,true(size(vok)));
        end
    end
else
    if isvector(v)
        vok = v(lok);
        vok = vok(:); % force column vectors
        doRepeated = false;
    else % v is a numeric matrix
        if isscalar(unique(sum(lok,2)))
            vok = v(lok);
            nCol = size(lok,2);
            vok = reshape(vok,[],nCol);
            % it is assumed the user provides repeated observations, every row
            % a subject and every column an observation
        else
            % x and y need to be reshaped into cells, every row a subject
            % and every cell a vector of replicates
            for r = flip(find(any(lok,2))).'
                vok{r,1} = v(r,lok(r,:)); %#ok<AGROW>
            end
        end
        doRepeated = true;
    end
end
end