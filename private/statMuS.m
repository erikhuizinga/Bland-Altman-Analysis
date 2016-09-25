function varargout = statMuS(varargin)
% mean-S statistics, S refers to either difference or ratio

%% determine statistic
SType = varargin{3};
switch SType
    case 'difference'
        % mean-difference
        SFun = @minus;
    case 'ratio'
        % mean-ratio
        SFun = @rdivide;
    case 'SD'
        % mean-standard deviation
        [x,y] = varargin{1:2};
        z = varargin{4};
        doConstantRegression = varargin{5};
        MSDType = varargin{6};
        varargout = statMuSD(x,y,z,doConstantRegression,MSDType);
        return
end
% SFun is the statistic function

%% input
% parse inputs
[x,y,~,n,z,t,doConstantRegression,assumeCTV] = varargin{:};
haveCell = iscell(x);

% vectorise x and y into X and Y and calculate the vectorised statistic
if haveCell
    X = [x{:}];
    Y = [y{:}];
else
    X = transpose(x);
    X = X(:);
    Y = transpose(y);
    Y = Y(:);
end
S = SFun(X,Y); % vectorised statistic
% X, Y and S now are row vectors of individual observations

% determine number of replicates
if haveCell
    m = cellfun(@numel,x);
else
    m = size(x,2)*ones(size(x,1),1);
end
% The number of replicates m is the same for x and y, because parseXY makes
% sure only pairs of observations are kept. m is an n×1 vector.

%% within subject means
% subject mean
if haveCell
    muXWithin = cellfun(@mean,x);
    muYWithin = cellfun(@mean,y);
else
    muXWithin = mean(x,2);
    muYWithin = mean(y,2);
end

% subject mean statistic, which is the statistic to plot as well
muSWithin = SFun(muXWithin,muYWithin); % SFun is the statistic function

%% calculate standard deviation of the statistic
if all(m==1)
    %% no repeated measurements, i.e. single mearuements
    % overall mean statistic, i.e. bias
    muS = mean(muSWithin);
    
    % variance of the statistic for single obsevations by the methods
    varS = var(muSWithin);
    varMuS = varS; % for single measurements
    
    % standard deviation of statistic for single obsevations by the methods
    sS = sqrt(varS);
    
    % within-subject variance is zero for single measurements
    varXWithin = 0;
    varYWithin = 0;
    
    % correction term factor is zero for single measurements
    corrM = 0;
else
    %% repeated measurements with equal/unequal numbers of replicates
    % prepare for ANOVA
    for sub = n:-1:1
        subjects{sub} = sub*ones(1,m(sub));
    end
    subjects = [subjects{:}]; % subject numbers, groups in ANOVA
    
    % total number of observation pairs
    N = numel(X); % equals numel(Y)
    
    % correction term factor
    corrM = 1-mean(1./m); % equals 1/n*sum(1./m), i.e. BA1999 notation
    
    % perform one-way ANOVA on X and Y with subject as groups
    for sub = n:-1:1 % loop over subjects
        % squared residual (SR) effect, per subject
        SRX(sub) = sum( ( X(subjects==sub)-muXWithin(sub) ).^2 );
        SRY(sub) = sum( ( Y(subjects==sub)-muYWithin(sub) ).^2 );
    end
    
    % SSR, sum of squared residual effects
    SSRX = sum(SRX);
    SSRY = sum(SRY);
    
    % MSR, mean squared residual effect
    dfR = N-n; % residual degrees of freedom
    MSRX = SSRX/dfR;
    MSRY = SSRY/dfR;
    
    % estimates of within-subject component variance
    varXWithin = MSRX;
    varYWithin = MSRY;
    
    if assumeCTV
        %% assume constant true value
        % mean statistic, i.e. bias
        muS = mean(SFun(muXWithin,muYWithin)); % BA1999§5.2
        
        % variance of the mean statistic for each subject
        varMuS = var(muSWithin);
        
        % variance of the statistic for single obsevations by the methods
        % BA1999 eq. 5.13, which reduces to eq. 5.2 for isscalar(unique(m))
        varS = varMuS + corrM*varXWithin + corrM*varYWithin;
        
        % standard deviation of the statistic for single obsevations by the methods
        sS = sqrt(varS);
    else
        %% assume variable true value
        % mean statistic, i.e. bias
        muS = mean(S); % according to BA1999§5.3 and BA2007§3
        
        % perform one-way ANOVA on statistic
        % SSS, sum of squared subjects effect
        SSSS = m'*( muSWithin-muS ).^2;
        
        % MSS, mean squared subjects effect
        MSSS = SSSS/(n-1);
        
        for sub = n:-1:1 % loop over subjects, groups in ANOVA
            % squared residual (SR) effect, per subject
            SRS(sub) = sum( ( S(subjects==sub)-muSWithin(sub) ).^2 );
        end
        
        % SSR, sum of squared residual effects
        SSRS = sum(SRS);
        
        % MSR, mean squared residual effect
        MSRS = SSRS/dfR;
        
        % estimates of within-subject component variance
        varSWithin = MSRS;
        
        % estimate of between-subject component variance
        n0 = (N^2-m'*m)/(n-1)/N; % see also Sahai 2005 p. 96
        varSBetween = (MSSS-MSRS)/n0;
        
        % estimate of total variance and standard deviation of the
        % statistic for single observations by the methods
        varS = varSBetween + varSWithin; % BA1999 p. 157 first eq.
        sS = sqrt(varS);
        
        % variance of the variance of components
        % The following two calculations are from the book Searle 2009. In
        % this function, the following variables refer to the following
        % book variables:
        % N: N
        % n: alpha
        % m: n
        % varSWithin: sigma_e^2
        % varSBetween: sigma_alpha^2
        varVarSWithin = 2*varSWithin^2/(N-n); % Searle 2006 p. 74 eq. 95 and Sahai 2005 p. 124 eq. 11.6.3
        
        % The book Searle 2006 Variance Components p. 75 eq. 102 gives:
        varVarSBetweenSearle2006 = ...
            2*N/( N^2 - m'*m ) ...
            * ( N*(N-1)*(n-1)/((N-n)*( N^2 - m'*m )) * varSWithin^2 + 2*varSWithin*varSBetween ...
            + ( N^2*(m'*m) + (m'*m)^2 - 2*N*sum(m.^3) )/( N*( N^2-m'*m ) ) * varSBetween^2 );
        % The formula of varVarSBetween above is written in the same order
        % as in Searle 2006 p. 75 eq. 102.
        
        %{
        % % The book Sahai 2005 gives three estimates from Crump and Searle:
        %
        % % Crump 1951 (from Sahai 2005 p. 124-5 eq. 11.6.4)
        % w = m*varSWithin./( varSWithin + m*varSBetween );
        % varVarSBetweenCrump1951 = ...
        %     2*varSWithin^2/n0^2 * ( 1/(n-1)^2 * (( 1/N*sum( m.^2./w.^2 ))^2 ...
        %     + sum( m.^2./w.^2 ) - 2/N*sum( m.^3./w.^3 )) + 1/(N-n) );
        % % Note: not the same (not even close) as
        % % varVarSBetweenSearle2006, which leads me to believe something's
        % % wrong in varVarSBetweenCrump1951 %FIXME
        %
        % % Searle 1956 (from Sahai 2005 p. 125 eq. 11.6.5)
        % tau = varSBetween/varSWithin;
        % varVarSBetweenSearle1956 = ...
        %     2*varSWithin^2/n0^2*( 1/( N*(n-1)^2 ) * ( tau^2/N * ( N^2*(m'*m) + (m'*m)^2 - 2*N*sum(m.^3) ) ...
        %     + 2*tau*( N^2 - m'*m )) + ( N-1 )/( (n-1)*(N-n)));
        % % Note: the same result as varVarSBetweenSearle2006 (logically)
        %
        % % Sahai 2005 also gives an alternative notation of the same Searle
        % % 1956 formula (Sahai 2005 p. 125 eq. 11.6.6)
        % varVarSBetweenSearle1956Alternative = ...
        %     2*varSWithin^2 * ( tau^2 * ( N^2*(m'*m) + (m'*m)^2 - 2*N*sum(m.^3) )/( N^2 - m'*m )^2 ...
        %     + 2*tau*N /( N^2 - m'*m ) + N^2*(N-1)*(n-1) /(( N^2 - m'*m )^2 * ( N-n )));
        % % Note: the same result as varVarSBetweenSearle2006 (logically)
        %
        % % Sahai 2005 p. 125 mentions another alternative:
        % % Rao 1997 p. 20 eq. 2.39
        % V = varSBetween + varSWithin./m;
        % varVarSBetweenRao1997 = ...
        %     2*N^2/( N^2 - m'*m )^2 ...
        %     *( sum( (m.^2).*( 1-2*m./N ).*V ) + ( (m.^2)'*V/N )^2 + (n-1)^2/(N-n)*varSWithin^2 );
        %}
        
        %%
        % variance of the global mean statistic, i.e. bias
        % From Sahai 2005 p. 102
        varMuS = m'*( varSWithin + m*varSBetween )/N^2;
        % sMuS = sqrt(varMuS);
    end
end

%% limits of agreement
loa = muS + z*sS*[-1 1];

%% subject mean to plot
muWithin = mean([muXWithin,muYWithin],2);

%% linear regression statistics
% linear regression of mean and S to plot with plotM
[polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
    linreg(muWithin,muSWithin,z,doConstantRegression);

%% confidence interval of the mean statistic (bias)
se2S = varS/n; % squared standard error of the bias
seS = sqrt(se2S); % bias standard error
eS = t*seS; % bias error
muSCI = muS + eS*[-1 1]; % bias confidence interval

%% confidence interval of the limits of agreement (loa)
if haveCell
    % unequal numbers of replicates
    varXW = cellfun(@var,x);
    varYW = cellfun(@var,y);
else
    % no or equal numbers of replicates
    varXW = var(x,[],2);
    varYW = var(y,[],2);
end
varVarXW = sum(2*varXW.^2./(m-1))/n^2; % BA1999 eq. 5.4, adjusted for unequal replicates
varVarYW = sum(2*varYW.^2./(m-1))/n^2;
% note these reduce to eq. 5.4 for equal replicates

% they are NaN in case all(mx==1)&all(my==1)
if isnan(varVarXW); varVarXW = 0; end
if isnan(varVarYW); varVarYW = 0; end

% BA1999 eq. 5.5 and 5.6 adjusted for unequal replicates
varCorrectionTerm = ...
    corrM^2 * varVarXW + ...
    corrM^2 * varVarYW;
% This is the variance of the correction term in eq. 5.3/5.13, adjusted
% for unequal replicates.

% BA1999 eq. 5.7
varVarMuS = 2*varMuS^2/(n-1);

% BA1999 eq. 5.8
varVarS = varVarMuS + varCorrectionTerm;

% BA1999 eq. 5.9
varSS = varVarS/4/varS;

% variance of the LOA
se2Loa = se2S + z^2*varSS; % BA1999 p. 153 eq. 5.10

% standard error of the LOA
seLoa = sqrt(se2Loa);

% confidence interval for the loa
eLoa = t*seLoa;
loaCI = [loa;loa] + eLoa*[-1 -1;1 1];

%% output
varargout = { ...
    muWithin,muSWithin, ...
    varXWithin,varYWithin, ...
    loaCI,loa, ...
    muS,muSCI, ...
    eLoa,eS, ...
    sS, ...
    polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa, ...
    m,X,Y ...
    };
end

function out = statMuSD(x,y,z,doConstantRegression,MSDType)
switch MSDType
    case 'difference'
        SFun = @minus;
    case 'ratio'
        SFun = @rdivide;
    case 'separate'
        SFun = @(x,y) x;
end

if iscell(x)
    % mean
    muX = cellfun(@mean,x);
    muY = cellfun(@mean,y);
    % statistic
    S = cellfun(SFun,x,y,'UniformOutput',false);
    % std of the statistic
    s = cellfun(@std,S);
else
    % mean
    muX = mean(x,2);
    muY = mean(y,2);
    % statistic
    S = SFun(x,y);
    % std of the statistic
    s = std(S,[],2);
end

mu = SFun(muX,muY);

[polyMSD,msePolyMSD,~,polyLLoaMSD,polyULoaMSD] = ...
    linreg(mu,s,z,doConstantRegression);

% output
out = {mu,s,polyMSD,msePolyMSD,polyLLoaMSD,polyULoaMSD};
end