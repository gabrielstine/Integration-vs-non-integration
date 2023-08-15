function [NLL, f] = GetDTB_FP_NLL(params, data, boundShape, stimStdev, fitFlag)
% Calculates the negative log-likelihood of data given parameters of
% a diffusion to bound model.  

% data takes the following form:

%   data = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (max nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (max nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)

% max nTrials is nTrials for the stim condition that has the most # of
% trials (usually 0% coherence). For the vectors RTs and choice, elements nTrials+1:maxnTrials
% should be NaNs. That way [data.RTs] or [data.choice] can be called as a
% matrix.

% boundShape is a string, specifying function to be used for calculating
% bounds. See expand_bounds.m

% fitFlag is 0 if you want to calculate things for the final output. Make
% 1 if fitting to make things go a little faster.

% Define params
k       = params(1) ;
b0      = params(2) ;
a       = params(3) ;
d       = params(4) ;
tNDmean = params(5) ;
tNDstd  = params(6) ;
offset  = params(7) ;
v       = params(8) ;

v = 0 ;

% v = 2 ;

% tNDmean = .4676 ;
% tNDstd  = .0867 ;

% if length(params) > 7
%     y0 = params(8) ;
% else
%     y0 = 0 ;
% end

dt = .0005 ;

mu = k*([data.coherence] - offset) ;
% var = 1+ (k*(stimStdev^2)) ;
var = 1 + (v*abs([data.coherence] - offset)) ;

RTs = [data.RTs] ;

tMax = 3 * max(max(RTs)) ;


% Create Bounds
t = 0:dt:tMax ;

Bup = expand_bounds(t,b0,a,d,boundShape)' ;
Blo = -Bup ;

% tMax = (max(RTs)) ;
% P = cell(1,length(mu)) ;
% tNDdist = P ;
% t = P ;

% parfor i = 1:length(mu)
%     
%     t{i} = 0:dt:tMax(i) ;
% 
%     Bup = expand_bounds(t{i},b0,a,d,boundShape)' ;
%     Blo = -Bup ;
%     
%     y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico
% 
%     y0 = zeros(size(y));
%     y0(findclose(y,0)) = 1;
%     y0 = y0/sum(y0);
%     
%     P{i} = dtb_fp_cc_vec(mu(i),t{i},Bup,Blo,y,y0,0,'var',var) ;
%     tNDdist{i} = dt * normpdf(t{i}, tNDmean, tNDstd)' ;
%     
% end

%%

% Solve FP to get decision time distributions
y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico

y0 = zeros(size(y));
y0(findclose(y,0)) = 1;
y0 = y0/sum(y0);

% P = cell(1,length(mu)) ;
% 
% parfor i = 1:length(mu)
%     P{i} = dtb_fp_cc_vec(mu(i),t,Bup,Blo,y,y0,0) ;
% end



P = dtb_fp_cc_vec(mu,t,Bup,Blo,y,y0,0,0,'var',var) ;
% toc
%%
% Create tND distribution
tNDdist = dt * normpdf(t, tNDmean, tNDstd)' ;


% Convolve tND distribution with DT distributions and calculate likelihood

% Preallocate
rtDistUp   = zeros(length(t)*2-1, length(mu)) ;
rtDistLo   = zeros(length(t)*2-1, length(mu)) ;
Likelihood = nan(size(RTs)) ;

% rtDistUp   = cell(1,length(mu)) ;
% rtDistLo   = cell(1,length(mu)) ;
% Likelihood = nan(size(RTs)) ;
    
% for i = 1:length(mu) 
%     
%     rtDistUp{i} = conv(tNDdist{i}, P{i}.up.pdf_t(1,:)') ;
%     rtDistLo{i} = conv(tNDdist{i}, P{i}.lo.pdf_t(1,:)') ;
%     
%     Likelihood(data(i).choice == 1, i) = rtDistUp{i}z(round(RTs(data(i).choice==1, i)/dt+1), i) ;
%     Likelihood(data(i).choice == 0, i) = rtDistLo(round(RTs(data(i).choice==0, i)/dt+1), i) ;
%     
% end

for i = 1:length(mu) 
    
    rtDistUp(:,i) = conv(tNDdist, P.up.pdf_t(i,:)') ;
    rtDistLo(:,i) = conv(tNDdist, P.lo.pdf_t(i,:)') ;
    
    Likelihood(data(i).choice == 1, i) = rtDistUp(round(RTs(data(i).choice==1, i)/dt+1), i) ;
    Likelihood(data(i).choice == 0, i) = rtDistLo(round(RTs(data(i).choice==0, i)/dt+1), i) ;
    
end

Likelihood(Likelihood < eps) = eps ;

NLL = -nansum(nansum(log(Likelihood))) ;


if ~fitFlag
    
%     tMax = 3 * max(max(RTs)) ;
% 
% 
% %     Create Bounds
%     t = 0:dt:tMax ;
% 
%     Bup = expand_bounds(t,b0,a,d,boundShape)' ;
%     Blo = -Bup ;
%     
%     y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico
%     y0 = zeros(size(y));
%     y0(findclose(y,0)) = 1;
%     y0 = y0/sum(y0);

    c4p = linspace(-1, 1, 50) ;
    mu4p = k*(c4p - offset) ;
    var2 = 1 + (v*abs(c4p - offset)) ;
    
    P2 = dtb_fp_cc_vec(mu4p,t,Bup,Blo,y,y0,0,0,'var',var2) ;
    
    muPrediction = k*([-1 1] - offset) ;
    Pprediction = dtb_fp_cc_vec(muPrediction,t,Bup,Blo,y,y0,0,0,'var',var) ;
    
    for i = 1:length(muPrediction)
        rtDistUpPred(:,i) = conv(tNDdist, Pprediction.up.pdf_t(i,:)') ;
        rtDistLoPred(:,i) = conv(tNDdist, Pprediction.lo.pdf_t(i,:)') ;
    end
    
    rtDistBothPrediction = rtDistUpPred + rtDistLoPred ;
    rtMeanPred = (Pprediction.up.p .* Pprediction.up.mean_t) + (Pprediction.lo.p .* Pprediction.lo.mean_t) + tNDmean ;
    
    rtDistBoth = rtDistUp + rtDistLo ;
    rtMean = (P2.up.p .* P2.up.mean_t) + (P2.lo.p .* P2.lo.mean_t) + tNDmean ;
    
    f.t = t ;
    f.c4p = c4p ;
    f.rtDistBoth = rtDistBoth(1:length(t), :) ;
    f.rtDistBoth100coh = rtDistBothPrediction ;
    f.rtMean100coh = rtMeanPred ;
    f.rtDistUp = rtDistUp(1:length(t), :) ;
    f.rtDistLo = rtDistLo(1:length(t), :) ;
    f.predChoice = P2.up.p ;
    f.predRTmean = rtMean ; 
    f.Bup = Bup ;
    f.Blo = Blo ;
    f.totArea = P.up.p + P.lo.p ; % Should be 1 for all coherences
    f.P = P ;
    
    if any(f.totArea) < .999
        warning('Total area under RT curve is less than 1. t vector is likely too short')
    end
    
else
    f = [];
end

end
    


    
