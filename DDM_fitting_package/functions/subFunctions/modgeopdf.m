function y = modgeopdf(n, p)

p  = p(:, 1:max(n)) ;  
pN = p(n) ;
q = 1-p ;
q2 = zeros(1, length(q)+1) ;
q2(1) = 1 ;
q2(2:end) = q ;
cumProd = cumprod(q2) ;


qProd = cumProd(n) ;

% qProd = zeros(1,length(n)) ;

% for i = 1:length(n)    
%     qProd(i) = prod(q(1:n(i)-1)) ;
% end


y = pN .* qProd ;

end