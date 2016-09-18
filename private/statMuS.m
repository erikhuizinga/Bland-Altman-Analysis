function [muXY,S,varXW,varYW,loaCI,loa,muS,muSCI,eLoa,eMuS,sS,polyMuS,msePolyMuS, ...
    sResPolyMuS,polyLLoa,polyULoa] = ...
    statMuS(x,y,SType,n,z,t,doConReg,repType)
% mean-S statistics
% S refers to either difference or ratio

switch SType
    case 'difference'
        fun = @minus;
    case 'ratio'
        assert(strcmp(repType,'none'),['Ratio statistics not ' ...
            'implemented for repeated measurements.'])
        fun = @rdivide;
end

% statistics without considering repeated measurements
muXY = mean([x(:), y(:)],2);
S = fun(x(:),y(:));

if repType
    % repeated measurements statistics
    
    % within subject means
    muXW = mean(x,2); % $\overline{X}$
    muYW = mean(y,2); % idem for Y
    
    % within subject mean statistic
    muSW = fun(muXW,muYW); % $\overline{d}$
    muS = mean(muSW); % mean statistic, i.e. bias
    varMuSW = var(muSW); % $s_\overline{d}^2$
    
    % within subject variances (zero for no replicates)
    varXW = var(x,[],2); % elements in $s_{xw}^2$ in BA1999
    varYW = var(y,[],2); % idem for y
    
    % number of replicates per subject
    if iscell(x) % and thus iscell(y)
        dbstack, keyboard %TODO verwijderen, testen
        mx = cellfun(@numel,x);
        my = cellfun(@numel,y);
    else % x and y numeric matrices or vectors
        mx = size(x,2)*ones(size(x,1),1); % $m_x$
        my = mx; % because this point can only be reached in equal number
        % of replicates
    end
    
    % weights of the correction term in BA1999 p. 155 eq. 5.13
    correctionTermXWeight = ( 1 - sum(1./mx)/n );
    correctionTermYWeight = ( 1 - sum(1./my)/n );
    % These weights reduce to those of BA1999 p. 151 eq. 5.3 for
    % isscalar(unique(mx)) & isscalar(unique(my)) and reduce to zero for
    % all(mx==1)&all(my==1), i.e. BA1999 section 2.
    
    % variance of the statistic
    % BA1999 p. 155 eq. 5.3 or 5.13
    varS = varMuSW + ...
        correctionTermXWeight * mean(varXW) + ...
        correctionTermYWeight * mean(varYW); % Var$(D)$
    sS = sqrt(varS);
    
    % loa of the statistic
    loa = muS + z*sS*[-1 1];
    
    % confidence interval of the bias
    varMuS = varS/n; % BA1999 p. 153 ‘variance of the mean difference $\overline{d}$’
    seMuS = sqrt(varMuS); % standard error of the bias
    eMuS = t*seMuS; % bias error
    muSCI = muS + eMuS*[-1 1];
    
    % % variance of the variance of S and of the standard deviation of S
    % varVarS = ... % BA1999 p. 152 eq. 5.8
    %     2*varMuSW^2/(n-1) + ...
    %     2*(mx-1)*varXW^2/n/mx^2 + ...
    %     2*(my-1)*varYW^2/n/my^2;
    
    % BA1999 p. 152 eq. 5.4, adjusted for unequal replicates
    varVarXW = sum(2*varXW.^2./(mx-1))/n/n; % $var(s_{xw}^2)$
    varVarYW = sum(2*varYW.^2./(my-1))/n/n; % idem for y
    % note these reduce to eq. 5.4 for equal replicates
    
    % these are nan in case all(mx==1)&all(my==1)
    if isnan(varVarXW); varVarXW = 0; end
    if isnan(varVarYW); varVarYW = 0; end
    
    % BA1999 p. 152 eq. 5.5 and 5.6 adjusted for unequal replicates
    varCorrectionTerm = ...
        correctionTermXWeight^2 * varVarXW + ...
        correctionTermYWeight^2 * varVarYW;
    % This is the variance of the correction term in eq. 5.3/5.13, adjusted
    % for unequal replicates.
    
    % BA1999 p. 152 eq. 5.7
    varVarMuS = 2*varMuSW^2/(n-1);
    
    % BA1999 p. 152 eq. 5.8
    varVarS = varVarMuS + varCorrectionTerm; % Var$(\hat{\sigma}_d^2)$
    
    % BA1999 p. 153 eq. 5.9
    varSS = varVarS/4/varS; % Var$(\hat{\sigma}_d)$
    
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