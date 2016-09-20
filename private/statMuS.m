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
[x,y,~,n,z,t,doConstantRegression] = varargin{:};
haveCell = iscell(x);

% vectorise x and y into X and Y
if haveCell
    X = [x{:}];
    Y = [y{:}];
else
    X = transpose(x);
    X = X(:);
    Y = transpose(y);
    Y = Y(:);
end

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
subjects = [subjects{:}]; % subject numbers, ‘groups’ in ANOVA

% subject mean
if haveCell
    muXW = cellfun(@mean,x);
    muYW = cellfun(@mean,y);
else
    muXW = mean(x,2);
    muYW = mean(y,2);
end

if all(m==1)
    varXWithin = 0;
    varYWithin = 0;
else
    % perform one-way ANOVA
    N = numel(X); % equals numel(Y)
    dfE = N-n; % error degrees of freedom
    for sub = n:-1:1 % loop over subjects, ‘groups’ in ANOVA
        % squared errors per subject
        SEX(sub) = sum( ( X(subjects==sub)-muXW(sub) ).^2 );
        SEY(sub) = sum( ( Y(subjects==sub)-muYW(sub) ).^2 );
    end
    SSEX = sum(SEX); % sum of squared errors
    SSEY = sum(SEY);
    varXWithin = SSEX/dfE; % estimate of within-subject variance for x, MSE
    varYWithin = SSEY/dfE; % estimate of within-subject variance for y, MSE
end

%% calculation of the statistic
% subject mean statistic, which is the statistic to plot as well
muSWithin = SFun(muXW,muYW); % SFun is the statistic function

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
mu = mean([muXW,muYW],2);

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
