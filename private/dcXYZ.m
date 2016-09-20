function str = dcXYZ(data,x,y,z)
str = {};
if ~isempty(data.xName), str = [str, [data.xName ': ' num2str(x)]]; end
if ~isempty(data.yName), str = [str, [data.yName ': ' num2str(y)]]; end
if ~isempty(z) && ~isempty(data.zName)
    str = [str, [data.zName ': ' num2str(z)]];
end
end
