function setOP(ax)
% make sure title is in visible figure area
op = ax.OuterPosition;
% not-so-nice attempt at detecting the subplot's maximum outerposition
maxPos = 1./(1:10); % possible positions
dPos = op(2)+op(4)-maxPos;
[~,iMin] = min(abs(dPos));
maxPos = maxPos(iMin);
ax.OuterPosition(4) = min(op(4)+op(2),maxPos) - op(2);