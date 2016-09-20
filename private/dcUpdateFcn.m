function str = dcUpdateFcn(~,e)
pos = e.Position;
x = pos(1);
y = pos(2);
if numel(pos) == 3;
    z = pos(3);
else
    z = [];
end
UserData = e.Target.UserData;
if iscell(UserData)
    data = UserData{1};
    m = UserData{2};
else
    data = UserData;
end
str = data.dcFun(data,x,y,z);
if ~isempty(data.prefix), str = [data.prefix, str]; end
if ~isempty(z) && exist('m','var')
    str = [str, ['m(i): ' num2str(m(z))]];
end
if ~isempty(data.suffix), str = [str, data.suffix]; end
end
