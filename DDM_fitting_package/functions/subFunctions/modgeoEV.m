function [ev, maxN] = modgeoEV(p) 

ev = zeros(1,size(p,1)) ;

for i = 1:size(p,1) 
    
    for j = 1000:1000:size(p,2)
        
        y = modgeopdf(j, p(i,:)) ;
        
        if y == 0
            maxN(i) = j ;
            break
        end
        
        if j == max(j)
            maxN(i) = max(j) ;
        end
        
    end
    
    ev(i) = sum(modgeopdf(1:maxN(i),p(i,:)) .* [1:maxN(i)]) ;

end


end