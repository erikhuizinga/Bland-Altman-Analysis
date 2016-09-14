function tf = validateXY(x)
% validate x,y input arguments of ba

if iscell(x)
    x = cellfun(@(x) x(:),x,'UniformOutput',false);
    x = vertcat(x{:});
end
tf = isnumeric(x);
end