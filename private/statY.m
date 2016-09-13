function [loaCI,loa,muY,muYCI,eLoaY,eMuY,sY,polyMuY,msePolyMuY, ...
    sResPolyMuY,polyLLoa,polyULoa] = statY(mu,y,n,z,t)
% Y statistics
% Y can be difference or ratio

muY = mean(y); % mean of difference or ratio

sY = std(y); % s_d in article
% varD = var(d); % s_d^2 in article

varMuY = sY^2/n; % variance of muD (SE^2) (p. 141)
seMuY = sqrt(varMuY); % standard error of muD (p. 142)

% varSD = sD^2/2/(n-1); % approximated variance of sD (article p. 141)
% Instead, use exact formula for unbiased estimated variance
% Source: http://stats.stackexchange.com/a/28567/80486
% sSD = sD * gamma((n-1)/2) / gamma(n/2) * ...
%     sqrt( (n-1)/2 - ( gamma(n/2) / gamma((n-1)/2) )^2 );
gammafrac = gamma((n-1)/2) / gamma(n/2);
% if ~isfinite(gammafrac) % true for large n
% % approximate using gamma(a+b)/gamma(a) ~ a^b
% % Source: https://en.wikipedia.org/w/index.php?title=Gamma_function&oldid=737220343#General
% % compare:
% % figure
% % n = 1:500;
% % g1 = gamma((n-1)/2)./gamma(n/2);
% % g2 = ((n-1)/2).^(-1/2);
% % g3 = sqrt(2./(n-1));
% % plot(n,g1,n,g2,n,g3,n,g1-g2,n,g1-g3,n,g2-g3)
% % legend g1 g2 g3 g1-g2 g1-g3 g2-g3
% gammafrac = sqrt(2/(n-1)); % same as: gammafrac = ((n-1)/2).^(-1/2);
% end
sSY = sY * gammafrac * sqrt( (n-1)/2 - gammafrac^-2 );
varSY = sSY^2; % unbiased estimate of variance of s_d

% limits of agreement statistics
loa = muY + z*sY*[-1 1]; % limits of agreement (LOA)

% confidence intervals (CI) for muD and loa
varLoa = varMuY + z^2*varSY; % article: bottom of p. 141
seLoa = sqrt(varLoa); % standard error of the LOA
eLoaY = t*seLoa; % LOA error
eMuY = t*seMuY; % muD error
muYCI = muY + eMuY*[-1 1];
loaCI = [loa;loa] + eLoaY*[-1 -1;1 1];
% loaCI is a 2x2 matrix. Every column contains the CI for a LOA: the first
% column corresponds to loa(1), the second to loa(2). The first and second
% row correspond to the lower and upper CI bound respectively.
% LOA = [loaCI(1,:);loa;loaCI(2,:)]; % optional 3x2 matrix form

% linear regression of mean and Y
[polyMuY,msePolyMuY,sResPolyMuY,polyLLoa,polyULoa] = baLinReg(mu,y,z);
end