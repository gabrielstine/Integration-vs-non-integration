function [NLL, f] = GetExtremaColBoundsNLL_sepRT(params, data, USfunc, stimStdev, fitFlag)

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
% tNDmean = .4676 ;
% tNDstd  = .0867 ;

dt = .0005 ;

tNDmean = tNDmean / dt ;
tNDstd  = tNDstd  / dt ;



k     = k*dt ;
sigma = dt + ([data.coherence]-offset)*v*dt ;
mu    = k * ([data.coherence] - offset) ;

if any(mu == 0)
    mu(mu==0) = 1e-30 ;
end

x  = 1:3 * max(max(round([data.RTs] ./ dt))) ;
Bs = expand_bounds(x*dt,b0,a,d,USfunc) ;

pRtemp = 1-normcdf(Bs,repmat(mu,length(Bs),1), repmat(sqrt(sigma),length(Bs),1))' ;
pLtemp = normcdf(-Bs,repmat(mu,length(Bs),1), repmat(sqrt(sigma),length(Bs),1))' ;

p  = pRtemp + pLtemp ;
pR = pRtemp ./ p ;



% Create tND distribution

tNDdist = normpdf(x, tNDmean, tNDstd) ;

Likelihood = nan(max([data.nTrials]), length(mu)) ;

pDT      = zeros(length(x), length(mu)) ;
pDTr     = zeros(length(x), length(mu)) ;
pDTl     = zeros(length(x), length(mu)) ;
rtDistR  = zeros(length(x)*2-1, length(mu)) ;
rtDistL  = zeros(length(x)*2-1, length(mu)) ;


for i = 1:length(mu)
    pDT(:,i)  = modgeopdf(x, p(i,:)) ;
    
%     pDTr(:,i) = (pR(i,:) .* pDT(:,i)') / sum(pR(i,:).*pDT(:,i)') ;
%     pDTl(:,i) = ((1-pR(i,:)) .* pDT(:,i)') / sum((1-pR(i,:)).*pDT(:,i)') ;
%     
    pDTr(:,i) = pR(i,:) .* pDT(:,i)' ;
    pDTl(:,i) = (1-pR(i,:)) .* pDT(:,i)' ;
    
    rtDistR(:,i) = conv(tNDdist, pDTr(:,i)) ;
    rtDistL(:,i) = conv(tNDdist, pDTl(:,i)) ;
    
    Likelihood(data(i).choice==1, i) = rtDistR(round(data(i).RTs(data(i).choice==1) ./ dt), i) ;
    Likelihood(data(i).choice==0, i) = rtDistL(round(data(i).RTs(data(i).choice==0) ./ dt), i) ;
end

Likelihood(Likelihood < eps) = eps ;

if all(all((~isnan([data.RTs])) == (~isnan(Likelihood))))
    NLL = -nansum(nansum(log(Likelihood))) ;
else
    NLL = 1e80 ;
end
    

if ~fitFlag

    
    ev = modgeoEV(p) ;
    pRight = nansum(pDTr) ;
    
    f.RTdistBoth  = rtDistR(1:length(x),:) + rtDistL(1:length(x),:) ;
    f.predDTr     = pDTr ;
    f.predDTl     = pDTl ;
    f.predChoice  = pRight ;
    f.predRTmean  = sum(x' .* f.RTdistBoth) * dt ;
    f.predRT_t    = x * dt ;
    f.predRTdistR = rtDistR(1:length(x), :) ;
    f.predRTdistL = rtDistL(1:length(x), :) ;
    f.bU = Bs ;
    f.bL = -Bs ;
else
    f = [];
end

end

