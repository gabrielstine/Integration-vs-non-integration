function [params, NLL, f] = FitExtremaVD(data, startValues, ntcRule, boundShape, varargin)

% Function to fit an extrema detection model to variable stimulus duration
% data. Minimizes the function GetVDExtremaColBoundsNLL. Similar
% inputs and outputs as FitAccumulationVD. 
%
% 
%   Parameter order is:
% 
%     k       = params(1) ;
%     b0      = params(2) ;
%     a       = params(3) ;
%     d       = params(4) ;
%     offset  = params(5) ;
%     guessBias = params(6) ;
%
%     ntcRule is a toggle for what to do when no threshold is crossed (i.e., no extremum is detected)
%     
%     if ntcRule == 1, use the sign of the last sample
%     if ntcRule == 2, flip a coin with bias equal to "guessBias"
% 
%         Guessing rule typically fits data much better
% 
%     Written by gms
%     Last updated 5/16/19

%%
% tic

k = 0 ;

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        k=varargin{i+1};
    end
    
end

options =   optimoptions(@fmincon, 'Display',   'iter',   'Maxiter',   1000,   'MaxFuneval',   2000,   'Algorithm',   'interior-point' ) ;

lBounds(1) = 20 ;
lBounds(2) = 0 ;
lBounds(3) = -50 ;
lBounds(4) = -50 ;
lBounds(5) =  -.5 ;
lBounds(6) =  0 ;

uBounds(1) = 200 ;
uBounds(2) = 1 ;
uBounds(3) = 50 ;
uBounds(4) = 50 ;
uBounds(5) = .5 ;
uBounds(6) = 1 ;


pl = [40 .04  .1 .01 -.1 .3]; 
ph = [150 .2  .6 .1   .1 .7];

% startValues = [12.1800    1.5200    2.2200    0.5100    0.2100    0.0900   0.0100] ;
% 
% boundShape = 'Logistic' ;

objFun  = @(params) GetVDExtremaColBoundsNLL(params, data, boundShape, ntcRule, 1, 'k',k) ;
params  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;
% params  = bads(objFun, startValues, lBounds, uBounds) ;

[NLL, f] = GetVDExtremaColBoundsNLL(params, data, boundShape, ntcRule, 0, 'k',k) ;
end