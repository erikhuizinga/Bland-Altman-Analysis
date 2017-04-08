function ZData = getZData(x, n, N, m)
if numel(x) == n
    ZData = 1 : n;
    
else  % numel(x) == N
    i = N;
    nn = n;
    mm = m;
    
    while i > 0
        ZData((i - mm(end) + 1) : i) = nn;
        i = i - mm(end);
        nn = nn - 1;
        mm(end) = [];
    end
end
end