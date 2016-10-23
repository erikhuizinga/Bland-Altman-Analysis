function str = dcUpdateFcn(~, e)
% Update data cursor string

% Get cursor coordinates
pos = e.Position;
x = pos(1);
y = pos(2);

if numel(pos) == 3
    z = pos(3);
else
    z = [];
end


% Get custom data
UserData = e.Target.UserData;
if iscell(UserData)
    data = UserData{1};
    m = UserData{2};
    
else
    data = UserData;
end


% Format new string
str = data.dcFun(data, x, y, z);

if ~isempty(data.prefix), str = [data.prefix, str]; end

if ~isempty(z) && exist('m', 'var')
    str = [str, ['m(i): ', num2str(m(z))]];
end

if ~isempty(data.suffix), str = [str, data.suffix]; end
end
