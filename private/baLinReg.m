function [polyXY,msePolyXY,sResPXY,polyLLoa,polyULoa] = ...
    baLinReg(x,y,z)
% simple linear regression statistics

%% stage 1
% simple linear regression of y on x
[polyXY,statsPolyXY] = polyfit(x,y,1); % article equation 3.1

% statistics
resPolyXY = y - polyval(polyXY,x); % residuals
ssePolyXY = sum(resPolyXY.^2); % or rss/ssr
msePolyXY = ssePolyXY/statsPolyXY.df;
% compare with the same calculation (CF toolbox required):
% [~,gof] = fit(x,y,'Poly1');
% gof.sse-ssePolyXY % equal within double precision
% gof.rmse-sqrt(msePolyXY) % idem
% gof.dfe-statsPolyXY.df % equal

sResPXY = std(resPolyXY); % standard deviation of the residuals

%% stage 2
% simple linear regression of the absolute residuals from stage 1 on x
R = abs(resPolyXY);
polyXRes = polyfit(x,R,1); % article equation 3.2

% limits of agreement (article equation 3.3)
polyLLoa = polyXY - z*sqrt(pi/2)*polyXRes; % lower LOA
polyULoa = polyXY + z*sqrt(pi/2)*polyXRes; % upper LOA

% If polyXRes(1) is not significantly different from zero, which is a
% matter clinical and statistical considerations, then the polynomial
% coefficients for the lower and upper LOA become
% polyLLoa = polyXY - z*sqrt(pi/2)*sResPXY
% polyULoa = polyXY + z*sqrt(pi/2)*sResPXY
end