function varargout = statMuS(varargin)
% mean-S statistics
% S refers to either difference or ratio

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
        varargout = statMuSD(x,y,z,doConstantRegression);
        return
end
% SFun is the statistic function

%% input
% parse inputs
[x,y,~,n,z,t,doConstantRegression,doCTV] = varargin{:};
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

% determine number of replicates
if haveCell
    m = cellfun(@numel,x);
else
    m = size(x,2)*ones(size(x,1),1);
end
% The number of replicates m is the same for x and y, because parseXY makes
% sure only pairs of observations are kept. m is an n×1 vector.

%% ANOVA
% prepare for one-way ANOVA
for sub = n:-1:1
    subjects{sub} = sub*ones(1,m(sub));
end
subjects = [subjects{:}]; % subject numbers, groups in ANOVA

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

% global mean statistic
muSGlobal = mean(S);

if all(m==1)
    varXWithin = 0;
    varYWithin = 0;
else
    % perform one-way ANOVA
    N = numel(X); % equals numel(Y)
    for sub = n:-1:1 % loop over subjects, groups in ANOVA
        % squared residual (SR) effect, per subject
        SRX(sub) = sum( ( X(subjects==sub)-muXWithin(sub) ).^2 );
        SRY(sub) = sum( ( Y(subjects==sub)-muYWithin(sub) ).^2 );
        SRD(sub) = sum( ( S(subjects==sub)-muSWithin(sub) ).^2 );
        % squared 
    end
    % SSS, sum of squared subjects effect
    SSSD = m.'*( muSWithin-muSGlobal ).^2;
    
    % MSS, mean squared subjects effect
    MSSS = SSSD/(n-1);
    
    % SSR, sum of squared residual effects
    SSRX = sum(SRX);
    SSRY = sum(SRY);
    SSRD = sum(SRD);
    
    % MSR, mean squared residual effect
    dfR = N-n; % residual degrees of freedom
    MSRX = SSRX/dfR;
    MSRY = SSRY/dfR;
    MSRS = SSRD/dfR;
    
    % estimates of within-subject component variance
    varXWithin = MSRX;
    varYWithin = MSRY;
    varSWithin = MSRS;
    
    % estimate of between-subject component variance
    divisor = (N*N-m'*m)/(n-1)/N;
    varSBetween = (MSSS-MSRS)/divisor;
    
    % estimate of total variance and standard deviation
    varTotalS = varSBetween + varSWithin;
    sTotalS = sqrt(varTotalS);
end

%% limits of agreement
% variance of the mean statistic
varMuS = var(muSWithin);

% correction term
corrM = 1-sum(1./m)/n; % equals mean(1./m), but this is BA1999's notation

% variance of the statistic
varS = varMuS + corrM*varXWithin + corrM*varYWithin;
% This is BA1999 eq. 5.13

% standard deviation of statistic for single obsevations by the methods
sS = sqrt(varS);

% overall mean statistic, i.e. bias
muS = mean(muSWithin);

% limits of agreement
loa = muS + z*sS*[-1 1];

%% calculation of the mean
% subject mean to plot
mu = mean([muXWithin,muYWithin],2);

%% linear regression statistics
% linear regression of mean and S to plot with plotM
[polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa] = ...
    baLinReg(mu,muSWithin,z,doConstantRegression);

%% confidence interval of the overall mean statistic (muS, i.e. bias)
se2S = varS/n; % squared standard error of the bias
seS = sqrt(se2S); % standard error
eS = t*seS; % bias error
muSCI = muS + eS*[-1 1]; % confidence interval

%% confidence interval of the limits of agreement (loa)
if haveCell
    varXW = cellfun(@(v) var(v,[],2),x);
    varYW = cellfun(@(v) var(v,[],2),y);
else
    varXW = var(x,[],2);
    varYW = var(y,[],2);
end
varVarXW = sum(2*varXW.^2./(m-1))/n/n; % BA1999 eq. 5.4, adjusted for unequal replicates
varVarYW = sum(2*varYW.^2./(m-1))/n/n;
% note these reduce to eq. 5.4 for equal replicates

% these are nan in case all(mx==1)&all(my==1)
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
    mu,muSWithin, ...
    varXWithin,varYWithin, ...
    loaCI,loa, ...
    muS,muSCI, ...
    eLoa,eS, ...
    sS, ...
    polyMuS,msePolyMuS,sResPolyMuS,polyLLoa,polyULoa, ...
    m,X,Y ...
    };
end

function out = statMuSD(x,y,z,doConstantRegression)
if iscell(x)
    % mean
    muX = cellfun(@mean,x);
    muY = cellfun(@mean,y);
    % std
    sX = cellfun(@std,x);
    sY = cellfun(@std,y);
else
    % mean
    muX = mean(x,2);
    muY = mean(y,2);
    % std
    sX = std(x,[],2);
    sY = std(y,[],2);
end

% regression of std on mean
[polyMSDX,msePolyMSDX,sResPolyMSDX,polyLLoaMSDX,polyULoaMSDX] = ...
    baLinReg(muX,sX,z,doConstantRegression);
[polyMSDY,msePolyMSDY,sResPolyMSDY,polyLLoaMSDY,polyULoaMSDY] = ...
    baLinReg(muY,sY,z,doConstantRegression);

% output
out = {muX,muY,sX,sY, ...
    polyMSDX,msePolyMSDX,sResPolyMSDX,polyLLoaMSDX,polyULoaMSDX, ...
    polyMSDY,msePolyMSDY,sResPolyMSDY,polyLLoaMSDY,polyULoaMSDY ...
    };
end
