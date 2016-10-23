function s2v(s)
% Convert structure fields to variables in caller workspace

%#ok<*NODEF>
fn = fieldnames(s);
for f = fn(:).'
    assignin('caller', f{:}, s.(f{:}))
end
end