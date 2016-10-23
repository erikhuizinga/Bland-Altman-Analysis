function setOP(ax)
% Make sure title is in visible figure area

% Get current outer position
op = ax.OuterPosition;

% Do a not-so-nice attempt at detecting the subplot's maximum outerposition
maxPos = 1 ./ (1 : 10);  % possible positions
dPos = op(2) + op(4) - maxPos;
[~, iMin] = min(abs(dPos));
maxPos = maxPos(iMin);

% Set new outer position
ax.OuterPosition(4) = min(op(4) + op(2), maxPos) - op(2);