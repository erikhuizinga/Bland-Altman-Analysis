function [x,y,doRepeated] = parseXY(x,y)

[x,repX] = parseV(x);
[y,repY] = parseV(y);

doRepeated = repX|repY;

% check number of elements
if isvector(x) && isvector(y)
    % x and y are numeric or cells
    if numel(x)~=numel(y)
        error 'Number of elements in x and y must be equal.'
    end
else % x and/or y matrices
    assert(size(x,1)==size(y,1),['x and y must have the same number of ' ...
        'rows (subjects) for repeated measures.'])
end
end

function [v,doRepeated] = parseV(v)
if iscell(v)
    assert(all(cellfun(@isvector,v)),['For cell input, every cell of ' ...
        inputname(1) ' must contain a vector.'])
    v = cellfun(@(v) v(:),v,'UniformOutput',false); % force column vectors
    nV = cellfun(@numel,v);
    if isscalar(unique(nV)), v = [v{:}]; end % convert to matrix
    doRepeated = true;
elseif isvector(v) % and not a cell
    v = v(:); % force column vectors
    doRepeated = false;
else % v is a numeric matrix
    doRepeated = true;
    % it is assumed the user provides repeated observations, every row a
    % subject and every column an observation
end
end