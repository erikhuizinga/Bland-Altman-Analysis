function [polyXY, msePolyXY, sResPXY, polyLLoa, polyULoa] ...
    = linreg(x, y, z, doConReg)
% Calculate simple linear regression statistics

%% Stage 1, calculate simple linear regression of y on x
% Regress y on x
[polyXY, statsPolyXY] = polyfit(x, y, 1);  % See BA1999 equation 3.1

% Calculate statistics
resPolyXY = y - polyval(polyXY, x);  % Calculate residuals (aka error)
ssePolyXY = sum(resPolyXY.^2);
msePolyXY = ssePolyXY / statsPolyXY.df;
% compare with the same calculation (CF toolbox required):
% [~, gof] = fit(x, y, 'Poly1');
% gof.sse - ssePolyXY  % equal within double precision
% gof.rmse - sqrt(msePolyXY)  % idem
% gof.dfe - statsPolyXY.df  % equal

% Calculate standard deviation of the residuals
sResPXY = std(resPolyXY);

%% stage 2, calculate simple linear regression of the absolute residuals from stage 1 on x
R = abs(resPolyXY);
polyXRes = polyfit(x, R, 1);  % See BA1999 equation 3.2

% Calculate limits of agreement, see BA1999 equation 3.3
if doConReg
    % If polyXRes(1) is not significantly different from zero, which is a
    % matter of clinical and statistical considerations, then the
    % polynomial coefficients for the lower and upper LOA become
    polyLLoa = polyXY - [0, z * sResPXY];
    polyULoa = polyXY + [0, z * sResPXY];
    
else
    polyLLoa = polyXY - z * sqrt(pi / 2) * polyXRes;  % lower LOA
    polyULoa = polyXY + z * sqrt(pi / 2) * polyXRes;  % upper LOA
end
end