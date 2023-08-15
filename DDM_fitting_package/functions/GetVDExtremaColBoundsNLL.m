function [NLL, f] = GetVDExtremaColBoundsNLL(params, data, USfunc, ntcRule, fitFlag, varargin)
%%
% Define params
if length(params) == 6
    k       = params(1) ;
    b0      = params(2) ;
    a       = params(3) ;
    d       = params(4) ;
    offset  = params(5) ;
    guessBias = params(6) ;

else
  
    k       = params(1) ;
    b0      = params(2) ;
    a       = params(3) ;
    d       = params(4) ;
    offset  = params(7) ;
end

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        if varargin{i+1} ~= 0
            k=varargin{i+1};
        end
    end

end

% k = 182.3656 ;

% guessBias = .5 ;

dt = .0005 ;

k     = k*dt ;
sigma = dt ;
mu    = k * ([data.coherence] - offset) ;

sdMat = [data.stimLength] ;
stimDurs = unique(sdMat(~isnan(sdMat))) / 1000 ;

if any(mu == 0)
    mu(mu==0) = 1e-30 ;
end


t = 1:3/dt;%max(stimDurs)/dt ;

thresh = expand_bounds(t*dt,b0,a,d,USfunc) ;


pR1 = 1-normcdf(thresh, mu, sqrt(sigma))' ;
pL1 = normcdf(-thresh, mu, sqrt(sigma))' ;

p  = pR1 + pL1 ;

pRtc  = pR1 ./ p ;

if ntcRule == 1 % last sample
    pR2 = (normcdf(thresh, mu, sqrt(sigma)) - normcdf(0, mu, sqrt(sigma)))' ;
    pL2 = (normcdf(0, mu, sqrt(sigma)) - normcdf(-thresh, mu, sqrt(sigma)))' ;
    pRntc = pR2 ./ (pR2 + pL2) ;
elseif ntcRule == 2 % guess
    pRntc = repmat(guessBias, length(t), length(mu))' ;
elseif ntcRule == 3 % Choose largest value so far (this will take real math)
    
%     for i = 1:length(mu)
%         for j = 1:length(stimDurs)
%             pRntc(i,j) = 
    
end


pTC = nan(length(mu),length(stimDurs)) ;
pDT = zeros(length(mu), length(t)) ;
pDTr = pDT ;
pDTl = pDT ;
nTrialsTDA = pTC ;
choiceData = pTC ;
realPR = pTC ;

% Get probability of crossing a threshold given stimulus duration and mu,
% and choice data

for i = 1:length(mu)
    
    pDT(i,:)  = modgeopdf(t, p(i,:)) ;
    pDTr(i,:) = pRtc(i,:) .* pDT(i,:) ;
    pDTl(i,:) = (1-pRtc(i,:)) .* pDT(i,:) ;
    
    for j = 1:length(stimDurs)
        
        nTrialsTDA(i,j) = length(data(i).choice(data(i).stimLength == stimDurs(j)*1000)) ;
        choiceData(i,j) = nansum(data(i).choice(data(i).stimLength == stimDurs(j)*1000)) ;
        
        pTC(i,j) = sum(pDT(i,1:round(stimDurs(j)/dt))) ;
        realPRtc(i,j) = sum(pDTr(i,1:round(stimDurs(j)/dt))) / (sum(pDTr(i,1:round(stimDurs(j)/dt))) + sum(pDTl(i,1:round(stimDurs(j)/dt)))) ;
        realPRntc(i,j) = pRntc(i, round(stimDurs(j)/dt)) ;
    end
    
end

predChoice = (pTC .* realPRtc ) + ((1-pTC).*realPRntc) ;
% predChoice(isnan(predChoice) = 

likelihood = binopdf(choiceData, nTrialsTDA, predChoice) ;
likelihood(likelihood < eps) = eps ;

% NLL = -nansum(nansum(log(likelihood))) ;
NLL = -sum(sum(log(likelihood))) ;
if isnan(NLL)
    NLL = 1e80 ;
end

if ~fitFlag
    
%     C4fit = linspace(.0001, .5, 25) ;
    
    f.stimDurs = stimDurs ;
    f.predChoice  = predChoice ;
    
    f.nTrialsTDA = nTrialsTDA ;
    f.predDTr = pDTr ;
    f.predDTl = pDTl ;
    f.bU = thresh ;
    f.bL = -thresh ;
    f.lapseRates = 1 - pTC ;
    
%     [~,simPredChoice,C4fit] = GetThresholds4ExtModPred(k, thresh, ntcRule, stimDurs) ;
%     [f.realThresholds, f.realThresholdsSEM] = GetThresholds(data) ;
    
    
    options = struct ;
    options.sigmoidName = 'weibull' ;
    options.expType = '2AFC' ;
    options.threshPC = .8 ;
    
    
    for i = 1:length(stimDurs)
        [f.realSens(:,i),~,stats] = glmfit([data.coherence], [choiceData(:,i) nTrialsTDA(:,i)], 'binomial', 'logit') ;
        f.predSens(:,i) = glmfit([data.coherence], [(round(predChoice(:,i)*30000))  repmat(30000,length(predChoice(:,i)),1)], 'binomial', 'logit') ;
        f.dataSensSE(i) = stats.se(2) ; 
        
%         sD = [C4fit' (round(simPredChoice(:,i)*1000)) repmat(1000,length(simPredChoice(:,i)),1)] ;
%         r  = psignifit(sD, options) ;
%         f.predThresholds(i) = exp(r.Fit(1)) ;
        
    end
    
    f.dataSens = f.realSens(2,:) ;
    f.predSens = f.predSens(2,:) ;
    
    
    

    f.choiceData = choiceData ./ nTrialsTDA ; 
else
    f = [];
end

end

