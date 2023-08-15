function [NLL, f] = GetDTB_VD_FP_NLL(params, data, boundShape, fitFlag, varargin)
% Calculates the negative log-likelihood of data given parameters of
% a drift diffusion model. Details of inputs and outputs can be found in
% FitAccumulationVD.


k       = params(1) ;
b0      = params(2) ;
a       = params(3) ;
d       = params(4) ;
cBias   = params(5) ;
dBias   = params(6) ;
% v       = params(7) ;

v = 0 ;

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        if varargin{i+1} ~= 0
            k=varargin{i+1};
        end
    end
    
    
    if isequal(varargin{i},'v')
        if varargin{i+1} ~= 0
            v=varargin{i+1};
        end
    end

end

% dBias = 0 ;

% k = 27.239  ;
% b0 = 1.0032 ;

dt = .0005 ;

C = [data.coherence] ;
mu    = k * (C - cBias) ;
var = 1 + (v*abs(C-cBias)) ;

sdMat = [data.stimLength] ;

stimDurs = unique(sdMat(~isnan(sdMat))) / 1000 ;

if any(mu == 0)
    mu(mu==0) = 1e-30 ;
end

t = 1:(max(stimDurs)/dt) + 10 ;

Bup = expand_bounds(t*dt,b0,a,d,boundShape) ;
Blo = -Bup ;

% Solve FP to get decision time distributions
y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico

y0 = normpdf(y,dBias,.01);
y0 = y0/sum(y0);

P = dtb_fp_cc_vec(mu,t*dt,Bup,Blo,y,y0,0,1,'var',var) ;

predChoice = P.up.p4VD(:,round(stimDurs/dt)) ;


for i = 1:length(mu)
    
    for j = 1:length(stimDurs)
        
        nTrialsTDA(i,j) = length(data(i).choice(data(i).stimLength == round(stimDurs(j)*1000))) ;
        choiceData(i,j) = nansum(data(i).choice(data(i).stimLength == round(stimDurs(j)*1000))) ;

    end
    
end

likelihood = binopdf(choiceData, nTrialsTDA, predChoice) ;
likelihood(likelihood < eps) = eps ;

NLL = -nansum(nansum(log(likelihood))) ;

choiceData2 = choiceData ./ nTrialsTDA ;



if ~fitFlag % Do extra calculations if you're no longer doing iterations
    
    f.params = params ;
    f.NLL    = NLL ;
    f.sCoh = C' ;
    f.stimDurs = stimDurs' ;
    f.nTrials  = nTrialsTDA ;
    f.predChoice = predChoice ;
    f.choiceData = choiceData2 ;
    f.bU = Bup ;
    f.bL = Blo ;
    
    t = 1:5/dt;

    Bup = expand_bounds(t*dt,b0,a,d,boundShape) ;
    Blo = -Bup ;


    % Solve FP to get decision time distributions
    y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico
    y0 = normpdf(y,dBias,.01);
    y0 = y0/sum(y0);
    P2 = dtb_fp_cc_vec(mu,t*dt,Bup,Blo,y,y0,0,1) ;
    
    
    C4t = linspace(.0001, .5, 25) ;
    mu4Thresh = k * (C4t - cBias) ;
    var4Thresh = 1 + (v * abs(C4t - cBias)) ;
    
    P4thresh = dtb_fp_cc_vec(mu4Thresh,t*dt,Bup,Blo,y,y0,0,1,'var',var4Thresh) ;
    
    pc4Thresh = P4thresh.up.p4VD(:, round(stimDurs/dt)) ;
    pc4Thresh(pc4Thresh>1) = 1 ;
% 
%     options = struct ;
%     options.sigmoidName = 'weibull' ;
%     options.expType = '2AFC' ;
%     options.threshPC = .632 ;
    
    for i = 1:length(stimDurs)

        f.predSens(:,i)           = glmfit(C4t, [round(pc4Thresh(:,i)*1000) repmat(1000,length(C4t),1)], 'binomial', 'logit') ;
        [f.dataSens(:,i),~,stats] = glmfit(C, [choiceData(:,i) nTrialsTDA(:,i)],'binomial','logit');
        
        predChoiceReg(:,i) = glmval(f.dataSens(:,i), C, 'logit') ;
        
        
      
        f.dataSensSE(i) = stats.se(2) ;
        
%         sD = [C4t' round(pc4Thresh(:,i)*1000) repmat(1000,length(C4t),1)] ;
%         r  = psignifit(sD, options) ;
%         f.predThresholds(i) = exp(r.Fit(1)) ;
        
        
    end
    
    LikeReg = binopdf(choiceData, nTrialsTDA, predChoiceReg) ;
    LikeReg(LikeReg<eps) = eps ;
    f.regNLL  = -sum(log(LikeReg(:))) ;
    
%     
    f.predSens = f.predSens(2,:) ;
    f.dataSens = f.dataSens(2,:) ;
    f.dt = dt ;
    f.t  = t ;
    f.P  = P2 ;
    f.P4thresh = pc4Thresh ;
    
    
%     if any(f.P.up.p + f.P.lo.p <.99)
%         warning('Probabilities do not sum to 1! Predicted decision times are likely off. Try increasing "t"')
%     end
    
else
    f = [];
end





