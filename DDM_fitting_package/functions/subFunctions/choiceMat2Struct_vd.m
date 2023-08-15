function d = choiceMat2Struct_vd(dataMat) 
%%
coh = dataMat(:,1);
usCoh = unique(coh) ;

for i = 1:length(usCoh)
    n2(i) = sum(coh == usCoh(i)) ;
end

maxnTrials = max(n2) ;

for i = 1:length(usCoh) 
    
    l = coh == usCoh(i) ;
    n = sum(l) ;
    
    d(i).nTrials = n ;
    d(i).coherence = usCoh(i) ;
    d(i).choice = nan(maxnTrials,1) ;
    d(i).stimLength = nan(maxnTrials,1) ;
    
    d(i).choice(1:n) = dataMat(l,2) ;
    d(i).stimLength(1:n) = dataMat(l,3) ;
    
end
    