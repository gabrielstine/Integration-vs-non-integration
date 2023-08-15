function d = choiceMat2Struct(dataMat) 
%% d = choiceMat2Struct(dataMat) 
% 
% Converts a data matrix w/ standard behavioral measurements into a structure
% to be used for fitting with Gabe's code. 
% 
% dataMat is an nTrials x 3 matrix with fields:
% 
% [signedCoherence    choice(0 for left, 1 for right)     RT(in seconds)]



coh = dataMat(:,1);
usCoh = unique(coh) ;

for i = 1:length(usCoh)
    n(i) = sum(coh == usCoh(i)) ;
end

maxnTrials = max(n) ;

for i = 1:length(usCoh) 
    
    l = coh == usCoh(i) ;
    n = sum(l) ;
    
    d(i).nTrials = n ;
    d(i).coherence = usCoh(i) ;
    d(i).choice = nan(maxnTrials,1) ;
    d(i).RTs = nan(maxnTrials,1) ;
    
    d(i).choice(1:n) = dataMat(l,2) ;
    d(i).RTs(1:n) = dataMat(l,3) ;
    
end
    