function [muXY,S,varXW,varYW,loaCI,loa,muS,muSCI,eLoa,eMuS,sS,polyMuS,msePolyMuS, ...
    sResPolyMuS,polyLLoa,polyULoa] = ...
    statMuS(x,y,SType,n,z,t,doConReg,doEqRep)
% mean-S statistics
% S refers to either difference or ratio

switch SType
    case 'difference'
        fun = @minus;
    case 'ratio'
        if doEqRep
            error(['Ratio statistics not implemented for repeated ' ...
                'measurements.'])
        end
        fun = @rdivide;
end

% statistics without considering repeated measurements
muXY = mean([x(:), y(:)],2);
S = fun(x(:),y(:));

if doEqRep
    % repeated measurements statistics
    
    % within subject means
    muXW = mean(x,2); % $\overline{X}$ in BA1999
    muYW = mean(y,2); % idem for Y
    
    % within subject mean statistic
    muSW = fun(muXW,muYW); % $\overline{D}$ in BA1999
    muS = mean(muSW); % mean statistic, i.e. bias
    varMuSW = var(muSW); % $s_\overline{d}^2$ in BA1999 p. 151
    
    % within subject variances
    varXW = mean(var(x,[],2)); % $s_{xw}^2$ in BA1999
    varYW = mean(var(y,[],2)); % idem for y
    
    % number of observations
    mx = size(x,2); % $m_x$ in BA1999
    my = size(y,2); % idem for y
    
    % variance of the statistic
    varS = varMuSW + (1-1/mx)*varXW + (1-1/my)*varYW; %BA1999 p. 151 eq. 5.3
    sS = sqrt(varS);
    
    % confidence interval of the bias
    varMuS = varS/n;
    seMuS = sqrt(varMuS);
    eMuS = t*seMuS; % muD error
    muSCI = muS + eMuS*[-1 1];
    
    % loa of the statistic
    loa = muS + z*sS*[-1 1];
    
    % variance of the variance of S and of the standard deviation of S
    varVarS = ... % BA1999 p. 152 eq. 5.8
        2*varMuSW^2/(n-1) + ...
        2*(mx-1)*varXW^2/n/mx^2 + ...
        2*(my-1)*varYW^2/n/my^2;
    varSS = varVarS/4/varS; % BA1999 p. 153 eq. 5.9
    
    % variance of the LOA
    varLoa = varMuS + z^2*varSS; % BA1999 p. 153 eq. 5.10
    
    % standard error of the LOA
    seLoa = sqrt(varLoa);
    
    % confidence interval for the loa
    eLoa = t*seLoa;
    loaCI = [loa;loa] + eLoa*[-1 -1;1 1];
    
    % linear regression of global mean and S
    [polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
        baLinReg(muXY,S,z,doConReg);
else % no repeated measurements, regular BAA LOA calculation
    % mean statistics
    muS = mean(S); % mean statistic, i.e. bias
    
    % standard deviation of the statistic
    sS = std(S); % $s_d$ in BA1999 p. 136
    
    varMuS = sS^2/n; % variance of muS, BA1999 p. 141
    seMuS = sqrt(varMuS); % standard error of muS, BA1999 p. 142
    
    % varSS = sS^2/2/(n-1); % approximated variance of sS, BA1999 p. 141
    % Instead, use exact formula for unbiased estimated variance
    % Source: http://stats.stackexchange.com/a/28567/80486
    % sSS = sS * gamma((n-1)/2) / gamma(n/2) * ...
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
    sSS = sS * gammafrac * sqrt( (n-1)/2 - gammafrac^-2 );
    varSS = sSS^2; % unbiased estimate of variance of s_d
    
    % limits of agreement statistics
    loa = muS + z*sS*[-1 1]; % limits of agreement (LOA)
    
    % confidence intervals (CI) for bias (muSCI) and LOA (loaCI)
    varLoa = varMuS + z^2*varSS; % BA1999 p. 141 (bottom)
    seLoa = sqrt(varLoa); % standard error of the LOA
    eLoa = t*seLoa; % LOA error
    eMuS = t*seMuS; % muS error
    muSCI = muS + eMuS*[-1 1];
    loaCI = [loa;loa] + eLoa*[-1 -1;1 1];
    % loaCI is a 2x2 matrix. Every column contains the CI for a LOA: the
    % first column corresponds to loa(1), the second to loa(2). The first
    % and second row correspond to the lower and upper CI bound
    % respectively.
    % LOA = [loaCI(1,:);loa;loaCI(2,:)]; % optional 3x2 matrix form
    
    % linear regression of mean and S
    [polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
        baLinReg(muXY,S,z,doConReg);
    
    % within-subject variance of x and y is always zero in regular BAA LOA
    varXW = 0;
    varYW = 0;
end
end