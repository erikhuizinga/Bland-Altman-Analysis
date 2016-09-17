function [xok,yok,n,repType] = prepRep(x,y,doRepeated)
% check and prepare x and y for repeated measurements

% check validity of x and y for calculations
lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);

if doRepeated % BAA for repeated measurements, determined by parseXY
    % distinguish the number of replicates
    if iscell(x) || iscell(y)
        % Unequal number of replicates. This is certain, because if
        % original input was a cell and contained equal number of
        % observations per subject, it was converted to a matrix by the
        % call to parseXY.m in ba.m.
        
        % check contents of every cell
        
        
        repType = 'unequal';
        n = numel(x);
    else % x and y are not cells, thus matrices
        % x and y must have the same dimensions:
        %  - The number of rows is the number of subjects, which must be
        %    equal for both x and y. This is checked by parseXY, so this is
        %    certain at this point.
        %  - The number of columns is the number of observations for all
        %    subjects. This must also be equal. However, there might be NaN
        %    or Inf elements in either x or y, which yields the possibility
        %    of BAA for repeated measurements with unequal number of
        %    replicates.
        
        % logical indices of elements to keep
        lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
        
        % check shape after selection of elements
        if isscalar(unique(sum(lok,2)))
            % x and y contain the same number of replicates per subject,
            % thus their
            % x and y remain matrices and the number of replicates
            % stays equal
            xok = x(lok);
            yok = y(lok);
            nCol = size(x,2);
            xok = reshape(xok,[],nCol);
            yok = reshape(yok,[],nCol);
            n = size(xok,1);
            
            % equal number of replicates
            repType = 'equal';
        else
            % x and y need to be reshaped into cells, every row a subject
            % and every cell a vector of replicates
            for r = flip(find(any(lok,2))).'
                xok{r,1} = x(r,lok(r,:));
                yok{r,1} = y(r,lok(r,:));
            end
            reptype = 'unequal';
            n = numel(xok);
        end
    end
else % no repeated measurements
    % keep only values that can be used in calculations
    lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
    n = nnz(lok);
    xok = x(lok);
    yok = y(lok);
    repType = 'none';
end
end