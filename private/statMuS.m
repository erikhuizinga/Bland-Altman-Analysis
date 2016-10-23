function varargout = statMuS(varargin)
% Calculate mean-S statistics, where S refers to either difference, ratio
% or standard deviation (SD)


%% Determine statistic
SType = varargin{3};
switch SType
    case 'difference'
        SFun = @minus;
        
    case 'ratio'
        SFun = @rdivide;
        
    case 'SD'
        [x, y] = varargin{1 : 2};
        z = varargin{4};
        doConstantRegression = varargin{5};
        MSDType = varargin{6};
        
        % Call function specific for the mean-SD graph
        varargout = statMuSD(x, y, z, doConstantRegression, MSDType);
        return
end


%% Parse inputs
[x, y, ~, n, z, t, doConstantRegression, assumeCTV] = varargin{:};
haveCell = iscell(x);

% Vectorise x and y into X and Y and calculate the vectorised statistic
if haveCell
    X = [x{:}];
    Y = [y{:}];
else
    X = transpose(x);
    X = X(:);
    Y = transpose(y);
    Y = Y(:);
end

% Calculate vectorised statistic
S = SFun(X, Y);  % SFun is the statistic function
% X, Y and S now are row vectors of individual observations

% Determine number of replicates
if haveCell
    m = cellfun(@numel, x);
else
    m = size(x, 2) * ones(size(x, 1), 1);
end
% The number of replicates m is the same for x and y, because parseXY makes
% sure only pairs of observations are kept. m is an n×1 vector.


%% Calculate within subject means
% Calculate subject mean per method
if haveCell
    muXWithin = cellfun(@mean, x);
    muYWithin = cellfun(@mean, y);
else
    muXWithin = mean(x, 2);
    muYWithin = mean(y, 2);
end

% Calculate subject mean statistic, which is the statistic to plot as well
muSWithin = SFun(muXWithin, muYWithin);


%% Calculate standard deviation of the statistic
if all(m == 1)  % No repeated measurements, i.e. single mearuements
    
    % Calculate overall mean statistic, i.e. bias
    muS = mean(muSWithin);
    
    % Calculate variance of the statistic for single observations
    varS = var(muSWithin);
    varMuS = varS;
    
    % Calculate standard deviation of statistic for single observation
    sS = sqrt(varS);
    
    % Set within-subject variance to zero for single observations
    varXWithin = 0;
    varYWithin = 0;
    
    % Set correction term factor to zero for single observations
    corrM = 0;
    
else  % BAA for repeated measurements with (un)equal numbers of replicates
    % Prepare for ANOVA
    for sub = n:-1:1
        subjects{sub} = sub * ones(1, m(sub));
    end
    subjects = [subjects{:}]; % subject numbers, groups in ANOVA
    
    % Calculate total number of observation pairs
    N = numel(X);  % Equals numel(Y)
    
    % Calculate correction term
    corrM = 1-mean(1./m);  % Equal to 1/n*sum(1./m), i.e. BA1999 notation
    
    % Perform one-way ANOVA on X and Y with subjects as groups
    for sub = n:-1:1  % Loop over subjects, the groups in ANOVA
        % Calculate squared errors (SE) per subject
        SEX(sub) = sum((X(subjects == sub) - muXWithin(sub)).^2);
        SEY(sub) = sum((Y(subjects == sub) - muYWithin(sub)).^2);
    end
    
    % Note that various sources use different conventions: sometimes the
    % errors are called residuals and the abbreviation is SR instead of SE.
    % The between subject terms are also called treatment, level, effect or
    % regression and the abbreviation therefore sometimes is SR, just like
    % with residuals, while it is something different. Here, E denotes
    % error (e.g. SSE is sum of squared error) and S denotes subject
    % effects (e.g. SSS is sum of squared subject effects).
    
    % Calculate SSE: sum of squared errors
    SSEX = sum(SEX);
    SSEY = sum(SEY);
    
    % Calculate MSE: mean squared error
    dfE = N-n;  % Set degrees of freedom of error
    MSEX = SSEX/dfE;
    MSEY = SSEY/dfE;
    
    % Estimate variance of within subjects component
    varXWithin = MSEX;
    varYWithin = MSEY;
    
    
    % Estimate variance of between subjects component
    if assumeCTV  % Assume constant true value
        % Calculate mean statistic, i.e. bias
        muS = mean(SFun(muXWithin, muYWithin));  % See BA1999§5.2
        
        % Calculate variance of the mean statistic for each subject
        varMuS = var(muSWithin);
        
        % Calculate variance of the statistic for single obsevations. See
        % BA1999 eq. 5.13, which reduces to eq. 5.2 for isscalar(unique(m))
        varS = varMuS + corrM*varXWithin + corrM*varYWithin;
        
        % Calculate standard deviation of the statistic for single
        % observations
        sS = sqrt(varS);
        
    else  % Assume variable true value
        % Calculate mean statistic, i.e. bias
        muS = mean(S);  % See BA1999§5.3 and BA2007§3
        
        % Perform one-way ANOVA on statistic
        % Calculate SSSS: sum of squared subjects effect of the statistic
        SSSS = m' * (muSWithin - muS).^2;
        
        % Calculate MSS: mean squared subjects effect
        MSSS = SSSS / (n - 1);
        
        % Loop over subjects, the groups in ANOVA
        for sub = n:-1:1
            % Calculate squared error (SE) per subject
            SES(sub) = sum((S(subjects == sub) - muSWithin(sub)).^2);
        end
        
        % Calculate SSE: sum of squared errors
        SSES = sum(SES);
        
        % MSE, mean squared errors
        MSES = SSES / dfE;
        
        % Estimate variance of within-subject component
        varSWithin = MSES;
        
        % Estimate of between-subject component variance
        n0 = (N^2 - m'*m) / (n - 1) / N;  % See Sahai 2005 p. 96
        varSBetween = (MSSS - MSES) / n0;
        
        % Estimate total variance and standard deviation of the statistic
        % for single observations
        varS = varSBetween + varSWithin;  % See BA1999 p. 157 first eq.
        sS = sqrt(varS);
        
        %% Estimate variance of the variance of components
        % The following two calculations are from the book Searle 2009. In
        % this function, the following variables refer to the following
        % book variables:
        % N: N
        % n: alpha
        % m: n
        % varSWithin: sigma_e^2
        % varSBetween: sigma_alpha^2
        varVarSWithin = 2 * varSWithin^2 / (N - n);  % See Searle 2006 p. 74 eq. 95 and Sahai 2005 p. 124 eq. 11.6.3
        
        % The book Searle 2006 Variance Components p. 75 eq. 102 gives:
        varVarSBetweenSearle2006 = ...
            2 * N / (N^2 - m'*m) ...
            * (N * (N - 1) * (n - 1) / ((N - n) * (N^2 - m'*m)) * varSWithin^2 + 2 * varSWithin * varSBetween ...
            + (N^2 * (m'*m) + (m'*m)^2 - 2 * N * sum(m.^3)) / (N * (N^2 - m'*m)) * varSBetween^2);
        % The formula above is written in the same order as in Searle 2006
        % p. 75 eq. 102.
        
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
        varMuS = m' * (varSWithin + m*varSBetween) / N^2;
        % sMuS = sqrt(varMuS);
    end
end


%% Calculate limits of agreement
loa = muS + z * sS * [-1, 1];


%% Calculate subject mean to plot
muWithin = mean([muXWithin, muYWithin], 2);


%% Calculate linear regression statistics of mean and S
[polyMuS, msePolyMuS, sResPolyMuS, polyLLoa, polyULoa] ...
    = linreg(muWithin, muSWithin, z, doConstantRegression);


%% Calculate confidence interval of the mean statistic (bias)
% Calculate squared standard error of the bias
se2S = varS / n;

% Calculate standard error of the bias
seS = sqrt(se2S);
eS = t * seS;  % bias error

% Calculate confidence interval of the bias
muSCI = muS + eS * [-1, 1];


%% Calculate confidence interval of the limits of agreement (loa)
% Calculate within-subject variance
if haveCell  % Calculate for unequal numbers of replicates
    varXW = cellfun(@var, x);
    varYW = cellfun(@var, y);
else  % Calculate for no or equal numbers of replicates
    varXW = var(x, [], 2);
    varYW = var(y, [], 2);
end

% Calculate variance of within-subject variance
varVarXW = sum(2 * varXW.^2 ./ (m - 1)) / n^2;  % See BA1999 eq. 5.4, adjusted for unequal replicates
varVarYW = sum(2 * varYW.^2 ./ (m - 1)) / n^2;
% Note these equations reduce to eq. 5.4 for equal replicates

% Handle NaN: they are NaN if all(mx == 1) & all(my == 1)
if isnan(varVarXW); varVarXW = 0; end
if isnan(varVarYW); varVarYW = 0; end

% BA1999 eq. 5.5 and 5.6 adjusted for unequal replicates
varCorrectionTerm = corrM^2 * varVarXW + ...
                    corrM^2 * varVarYW;
% This is the variance of the correction term in eq. 5.3/5.13, adjusted
% for unequal replicates.

% See BA1999 eq. 5.7
varVarMuS = 2 * varMuS^2 / (n - 1);

% See BA1999 eq. 5.8
varVarS = varVarMuS + varCorrectionTerm;

% See BA1999 eq. 5.9
varSS = varVarS / 4 / varS;

% Calculate variance of the LOA
se2Loa = se2S + z^2 * varSS;  % See BA1999 p. 153 eq. 5.10

% Calculate standard error of the LOA
seLoa = sqrt(se2Loa);

% Calculate confidence interval of the loa
eLoa = t * seLoa;
loaCI = [loa; loa] + eLoa * [-1, -1 ; 1, 1];


%% Set output
varargout = {muWithin, muSWithin, varXWithin, varYWithin, loaCI, loa, ...
             muS, muSCI, eLoa, eS, sS, polyMuS, msePolyMuS, ...
             sResPolyMuS, polyLLoa, polyULoa, m, X, Y};
end


function out = statMuSD(x, y, z, doConstantRegression, MSDType)
switch MSDType
    case 'difference'
        SFun = @minus;
        
    case 'ratio'
        SFun = @rdivide;
        
    case 'separate'
        SFun = @(x, y) x;
end

if iscell(x)
    % Calculate mean
    muX = cellfun(@mean, x);
    muY = cellfun(@mean, y);
    
    % Calculate statistic
    S = cellfun(SFun, x, y, 'UniformOutput', false);
    
    % Calculate std of the statistic
    s = cellfun(@std, S);
    
else
    % Calculate mean
    muX = mean(x, 2);
    muY = mean(y, 2);
    
    % Calculate statistic
    S = SFun(x, y);
    
    % Calculate std of the statistic
    s = std(S, [], 2);
end

mu = SFun(muX, muY);

[polyMSD, msePolyMSD, ~, polyLLoaMSD, polyULoaMSD] ...
    = linreg(mu, s, z, doConstantRegression);

% Set output
out = {mu, s, polyMSD, msePolyMSD, polyLLoaMSD, polyULoaMSD};
end