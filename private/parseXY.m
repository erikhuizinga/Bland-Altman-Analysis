function [x,y,doRepeated] = parseXY(x,y)

[x,repX] = parseV(x);
[y,repY] = parseV(y);

doRepeated = repX|repY;

% check number of elements
if isvector(x) && isvector(y)
    % x and y are numeric or cells
    assert(numel(x)==numel(y), ...
        'Number of elements in x and y must be equal.')
else % x and/or y matrices
    % validate number of rows (subjects)
    assert(size(x,1)==size(y,1),['x and y must have the same number of ' ...
        'rows (subjects) for repeated measures.'])
    % check number of columns (replicates)
        assert(size(x,2)==size(y,2),['x and y must have the same ' ...
            'number of columns (replicates) per subject.'])
end
end

function [v,doRepeated] = parseV(v)
if iscell(v)
    doRepeated = true;
    assert(all(cellfun(@isvector,v)),['For cell input, every cell of ' ...
        inputname(1) ' must contain a vector.'])
    v = cellfun(@(v) v(:).',v,'UniformOutput',false); % force row vectors
    v = v(:); % force column cell
    % now v is a column vector cell, every element a subject and every
    % cell's contents a row vector of replicates
    nV = cellfun(@numel,v);
    if isscalar(unique(nV))
        % all cells contain the same number of replicates, thus they fit in
        % a matrix
        v = vertcat(v{:}); % convert to matrix
        % now v is a numeric matrix with rows corresponding to subjects and
        % columns corresponding to replicates
        if isvector(v)
            % a single measurement per cell element, thus no replicates,
            % reparse
            [v,doRepeated] = parseV(v);
        end
    end
elseif isvector(v)
    v = v(:); % force column vectors
    doRepeated = false;
else % v is a numeric matrix
    doRepeated = true;
    % it is assumed the user provides repeated observations, every row a
    % subject and every column an observation
end
end