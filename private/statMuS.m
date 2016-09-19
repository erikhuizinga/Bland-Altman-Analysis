function varargout = statMuS(varargin)
% mean-S statistics
% S refers to either difference or ratio

SType = varargin{3};
switch SType
    case 'difference'
        % mean-difference
        SFun = @minus;
    case 'ratio'
        % mean-ratio
        repType = varargin{8};
        assert(strcmp(repType,'none'),['Ratio statistics not ' ...
            'implemented for repeated measurements.'])
        SFun = @rdivide;
    case 'SD'
        % mean-standard deviation
        [x,y] = varargin{1:2};
        varargout = statMuSD(x,y);
        return
end
% SFun is the statistic function

% parse inputs
[x,y,~,n,z,t,doConstantRegression,~] = varargin{:};

if iscell(x) % and y is a cell too
    % first determine size of each cell    
    % number of replicates per subject
    % mx = cellfun(@numel,x);
    % my = mx; % this is enforced by parseXY
    m = cellfun(@numel,x);
    
    % prepare for one-way ANOVA
    for i = n:-1:1
        subjects{i} = i*ones(1,m(i));
    end
    subjects = [subjects{:}]; % subject numbers
    X = [x{:}]; % vectorise data
    Y = [y{:}]; % idem
    muXW = cellfun(@mean,x);
    muYW = cellfun(@mean,y);
    
    % perform one-way ANOVA
    N = numel(X);
    dfE = N-n;
    for sub = n:-1:1
        SSEX(sub) = sum( ( X(subjects==sub)-muXW(sub) ).^2 );
        SSEY(sub) = sum( ( Y(subjects==sub)-muYW(sub) ).^2 );
    end
    SSEX = sum(SSEX);
    SSEY = sum(SSEY);
    varXWithin = SSEX/dfE; % estimate of within-subject variance for x, MSE
    varYWithin = SSEY/dfE; % estimate of within-subject variance for y, MSE
    
    % correction term 
    corrM = sum(1./m)/n; % equal to mean(1./m), but this is the notation in BA1999
    
    % subject mean statistic
    muSW = SFun(muXW,muYW);
    
    % variance of the mean statistic
    varMuS = var(muSW);
    
    % variance of the statistic
    varS = varMuS + (1-corrM)*varXWithin + (1-corrM)*varYWithin;
    % This is BA1999 eq. 5.13
    
    % standard deviation of statistic for single obsevations by the methods
    sS = sqrt(varS);
    
    % overall mean statistic, i.e. bias
    muS = mean(muSW);
    
    % limits of agreement
    eMuS = z*sS;
    loa = muS + eMuS*[-1 1];
    
    % statistics without consideration of repeated measurements
    muXY = mean([X;Y]);
    S = SFun(X,Y);
    
    % subject mean to plot
    mu = mean([muXW,muYW],2);
    
    % linear regression of mean and S to plot with plotM
    [polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
        baLinReg(mu,muSW,z,doConstantRegression);
else
    % within subject means
    muXW = mean(x,2); % $\overline{X}$
    muYW = mean(y,2); % idem for Y
    
    % within subject mean statistic
    muSW = SFun(muXW,muYW); % $\overline{d}$
    
    % overall mean statistic
    muS = mean(muSW); % mean statistic, i.e. bias
    varMuSW = var(muSW); % $s_\overline{d}^2$
    
    % within subject variances (zero for no replicates)
    varXW = var(x,[],2); % elements in $s_{xw}^2$ in BA1999
    varYW = var(y,[],2); % idem for y
    
    mx = size(x,2)*ones(size(x,1),1); % $m_x$
    my = mx; % because this point can only be reached in equal number
    % of replicates
    
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
    
    % statistics without considering repeated measurements
    muXY = mean([x(:), y(:)],2);
    S = SFun(x(:),y(:));
    
    % linear regression of global mean and S
    [polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
        baLinReg(muXY,S,z,doConstantRegression);
end

varargout = { ...
    muXY,S, ...
    varXW,varYW, ...
    loaCI,loa, ...
    muS,muSCI, ...
    eLoa,eMuS, ...
    sS, ...
    polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa ...
    };

%% archive:
% early version of exact calculation of sSS and varSS
% % varSS = sS^2/2/(n-1); % approximated variance of sS, BA1999 p. 141
% % Instead, use exact formula for unbiased estimated variance
% % Source: http://stats.stackexchange.com/a/28567/80486
% % sSS = sS * gamma((n-1)/2) / gamma(n/2) * ...
% %     sqrt( (n-1)/2 - ( gamma(n/2) / gamma((n-1)/2) )^2 );
% gammafrac = gamma((n-1)/2) / gamma(n/2);
% % if ~isfinite(gammafrac) % true for large n
% % % approximate using gamma(a+b)/gamma(a) ~ a^b
% % % Source: https://en.wikipedia.org/w/index.php?title=Gamma_function&oldid=737220343#General
% % % compare:
% % % figure
% % % n = 1:500;
% % % g1 = gamma((n-1)/2)./gamma(n/2);
% % % g2 = ((n-1)/2).^(-1/2);
% % % g3 = sqrt(2./(n-1));
% % % plot(n,g1,n,g2,n,g3,n,g1-g2,n,g1-g3,n,g2-g3)
% % % legend g1 g2 g3 g1-g2 g1-g3 g2-g3
% % gammafrac = sqrt(2/(n-1)); % same as: gammafrac = ((n-1)/2).^(-1/2);
% % end
% sSS = sS * gammafrac * sqrt( (n-1)/2 - gammafrac^-2 );
% varSS = sSS^2; % unbiased estimate of variance of s_d
end

function out = statMuSD(x,y)
% mean
muX = mean(x,2);
muY = mean(y,2);
% std
sX = std(x,[],2);
sY = std(y,[],2);
% output
out = {muX,muY,sX,sY};
end