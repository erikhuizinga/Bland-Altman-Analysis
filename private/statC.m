function [pRhoXY, rhoXY, polyXY, msePXY] = statC(xok, yok, z, doConReg)
% Calculate correlation statistics

% Calculate correlation
[rhoXY, pRhoXY] = corrcoef(xok, yok);
rhoXY = rhoXY(1, 2);
pRhoXY = pRhoXY(1, 2);

% Calculate linear regression statistics
[polyXY, msePXY] = linreg(xok, yok, z, doConReg);
end