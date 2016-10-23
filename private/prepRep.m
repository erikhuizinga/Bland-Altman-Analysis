function [x, y, n] = prepRep(x, y, doRepeated)
% Check and prepare x and y for repeated measurements

if doRepeated  % BAA for repeated measurements, determined by parseXY
    % Distinguish between the options of number of replicates
    if iscell(x) || iscell(y)
        % Unequal number of replicates. This is certain, because if
        % original input was a cell and contained equal number of
        % observations per subject, it was converted to a matrix by the
        % call to parseXY.m in ba.m.
        n = numel(x);
        
    else  % x and y are not cells, thus matrices
        % Equal number of replicates
        % x and y must have the same dimensions:
        %  - The number of rows is the number of subjects, which must be
        %    equal for both x and y. This is checked by parseXY, so this is
        %    certain at this point
        %  - The number of columns is the number of observations for all
        %    subjects. This must also be equal and is checked by parseX.
        n = size(x, 1);
    end
    
else  % No repeated measurements
    n = numel(x);
end
end