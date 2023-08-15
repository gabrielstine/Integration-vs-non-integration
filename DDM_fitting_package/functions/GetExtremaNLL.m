function [NLL, f] = GetExtremaNLL(params, d, stimStdev, useChoiceFlag, varargin)

% Calculates the negative log-likelihood of data in "d" given parameters of
% an extrema aka "probability summation" model. 

% d takes the following form:

%   d = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
%
%
%
%   useChoiceFlag takes the value 0, 1, or 3 (don't ask)
%       if 0, calculate NLL using only RTs
%       if 1, calculate NLL using both RTs and choice data
%       if 3, calculate NLL using only choice data. Remember that you
%       cannot predict RTs with fits to choice data.



%   Written by GMS. Updated 7/5/17

%%

dT = 0.0005 ;
tNDmax_right = 5 ;
tNDmax_left  = 5 ;

% Define params
k         = params(1) ; % Kappa
A         = params(2) ; % Bound
tND_right = params(3) ; % non-decision time, righward
tND_left  = params(4) ; % non-decision time, righward
offset    = params(5) ; % Bias
% lRate     = params(6) ; % lapse rate
v         = params(6) ; % Variance parameter (0 if no scaling)

v = 0 ;

lRate = 0 ;

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        if varargin{i+1} ~= 0
            k=varargin{i+1};
        end
    end
    

    if isequal(varargin{i},'tNDmax_right')
        if varargin{i+1} ~= 0
            tNDmax_right=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tND_right')
        if varargin{i+1} ~= 0
            tND_right=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tNDmax_left')
        if varargin{i+1} ~= 0
            tNDmax_left=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tND_left')
        if varargin{i+1} ~= 0
            tND_left=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tNDmax')
        if varargin{i+1} ~= 0
            tNDmax_right=varargin{i+1};
            tNDmax_left=varargin{i+1};
        end
    end
    
end

tND_left = tND_right ;

% k = 122.38
% tND = 0.5027
k      = k*dT ;
v      = v*dT ;
% sigma  = dT + (k*stimStdev)^2 ; %+ (k*stimStdev)^2 ;
% sigma = dT;

% A = A*dT ;



% Fit glm to get bias term. Bias term will be used to select trials where
% choice matched the sign of Mu.

if length([d.coherence]) > 2

    % Create predictor matrix
    for i = 1:length(d)
        temp{i}  = repmat(d(i).coherence, length(d(i).RTs), 1) ;
    end
    
    glmX = cat(1, temp{:}) ;                            % Vector containing the coherence on each trial
    glmY = reshape([d.choice],numel([d.choice]), 1) ;   % Vector containing responses on each trial
    
    % Fit glm
    Bs   = glmfit(glmX, glmY, 'binomial', 'logit') ;
    bias = Bs(1) / Bs(2) ;
    
else
    bias = 0 ;
end
    
    
% Get reaction times for trials where choice matches Mu (i.e. rightward
% choices when C > bias and leftward choices when C < bias). Calc mean and
% variance.
for i = 1:length(d) 
    
    if d(i).coherence >= -bias
        goodRTs{i} = d(i).RTs(d(i).choice == 1) ;
    else
        goodRTs{i} = d(i).RTs(d(i).choice == 0) ;
    end
    
    meanRTs(i)   = nanmean(goodRTs{i}) ;
    predRTsem(i) = sqrt(nanvar(goodRTs{i}) ./ length(goodRTs{i})) ; 
    
end
    
% 
% if useChoiceFlag == 0 || useChoiceFlag == 3
%     offset = -bias ;
% end

mu  = k * ([d.coherence] - offset) ;
sigma  = dT + (v * abs([d.coherence] - offset)) ;



% Make sure tND can't be above some maximum
if tND_right > tNDmax_right
    tND_right = tNDmax_right ;
end

if tND_left > tNDmax_left
    tND_left = tNDmax_left ;
end




% Change zero coherence to very small to avoid using limits
if any(mu == 0)
    mu(mu==0) = 1e-30 ;
end

%%

% Calculate model predictions

    % Get RT prediction
    pR = 1 - normcdf(A, mu, sqrt(sigma)) ;     % Area under the curve beyond A
    pL = normcdf(-A, mu, sqrt(sigma)) ;        % Area under the curve beyond -A
    p  = pR + pL ;                             % Probability of momentary evidence crossing +-A
    
    predDT     = (1 ./ p) * dT ;               % Expected value of geometric dist. is 1/p
%     predDTvar  = ((1-p) ./ p.^2) * dT^2 ;      % Equation for variance of a geometric distribution

    leftMu = mu < 0 ;
    predRTmean(leftMu)  = predDT(leftMu)  + tND_left ;
    predRTmean(~leftMu) = predDT(~leftMu) + tND_right ;
%     predRTvar  = predDTvar + varTND ; % Variances add
%     predRTsem  = sqrt(predRTvar ./ [d.nTrials]) ;


    % Get accuracy prediction
    predChoice = lRate + (1 - 2 * lRate) .* (pR ./ p) ;
    
    
%%
% Calculate log-likelihood of observed mean RTs, given predicted mean RTs
% and predicted sem of RTs 
likelihoodRT = normpdf(meanRTs, predRTmean, predRTsem) ;
likelihoodRT(likelihoodRT == 0) = eps ;

rtNLL = -sum(log(likelihoodRT)) ;

% Get likelihood for choice
likelihoodChoice = binopdf(nansum([d.choice]), [d.nTrials], predChoice) ;
likelihoodChoice(likelihoodChoice == 0) = eps ;

choiceNLL = -sum(log(likelihoodChoice)) ;

% Calculate total NLL
if useChoiceFlag == 1
    NLL = rtNLL + choiceNLL ;
elseif useChoiceFlag == 0
    NLL = rtNLL ;
elseif useChoiceFlag == 3
    NLL = choiceNLL ;
end
    
if NLL == inf || isnan(NLL)
    NLL = 1e80 ;
end

% c4p = linspace(d(1).coherence, d(end).coherence, 100) ;
c4p = -1.05:.005:1.05 ;
mu4p = k*(c4p - offset) ;
v4p  = dT + (v*abs(c4p - offset)) ;

% mu4p = mu4p ./ v4p ;

lMu = mu4p < 0 ;

% Get RT prediction
pR = 1 - normcdf(A, mu4p, sqrt(v4p)) ;     % Area under the curve beyond A
pL = normcdf(-A, mu4p, sqrt(v4p)) ;        % Area under the curve beyond -A
p  = pR + pL ;                             % Probability of momentary evidence crossing +-A

predDT     = (1 ./ p) * dT ;               % Expected value of geometric dist. is 1/p
predRTmean(lMu)  = predDT(lMu)  + tND_left ;
predRTmean(~lMu) = predDT(~lMu) + tND_right ;

predChoice = lRate + (1 - 2 * lRate) .* (pR ./ p) ;


f.c4p = c4p ;
f.predChoice = predChoice ;
temp = nansum([d.choice]) ./ [d.nTrials] ;
f.choiceData = nansum([d.choice]) ./ [d.nTrials] ;
f.meanRTs = meanRTs ;
f.semRTs  = predRTsem ;
f.choiceSE   = sqrt( (temp.*(1-temp)) ./ [d.nTrials]) ;
f.predRTmean = predRTmean ;
% f.predRTsem  = predRTsem ;
f.bias = -bias ;
f.sigma = sigma ;
f.rtNLL = rtNLL ;




end