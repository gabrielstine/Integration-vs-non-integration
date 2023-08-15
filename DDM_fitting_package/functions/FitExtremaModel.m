function [params NLL f] = FitExtremaModel( d, startValues, useChoiceFlag, varargin)

% Function to fit RT and/or choice data using an extrema (aka "probability
% summation") model. Minimizes the function GetExtremaNLL.
% 
% See help for FitAccumulationModel, which has the same exact parameterization
% 
% Written by gms
% Last updated 5/16/19

k = 0 ;
tND_right = 0 ;
tND_left = 0 ;
tNDmax_right = 0 ;
tNDmax_left  = 0 ;
tNDmax = 0 ;

options =   optimset( 'Display',   'off',   'Maxiter',   1000,   'MaxFuneval',   4000,   'Algorithm',   'interior-point' ) ;

lBounds(1) = 5 ;
lBounds(2) = 0 ;
lBounds(3) = 0 ;
lBounds(4) = 0 ;
lBounds(5) = -.1 ;
lBounds(6) = 0 ;


uBounds(1) = 200 ;
uBounds(2) = .11 ;
uBounds(3) = 1 ;
uBounds(4) = 1 ;
uBounds(5) = .1 ;
uBounds(6) = 200 ;


for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        k=varargin{i+1};
    end
    
    if isequal(varargin{i},'tND_right')
        tND_right=varargin{i+1};
    end
    
    if isequal(varargin{i},'tNDmax_right')
        tNDmax_right=varargin{i+1};
    end

    if isequal(varargin{i},'tND_left')
        tND_left=varargin{i+1};
    end
    if isequal(varargin{i},'tNDmax_right')
        tNDmax_right=varargin{i+1};
    end

    if isequal(varargin{i},'tND_left')
        tND_left=varargin{i+1};
    end
    
    if isequal(varargin{i},'tNDmax_left')
        tNDmax_left=varargin{i+1};
    end    

    if isequal(varargin{i},'tNDmax')
        tNDmax=varargin{i+1};
    end 
    
end




objFun  = @(params) GetExtremaNLL(params, d, 0, useChoiceFlag, ...
    'k',k,'tND_right',tND_right,'tNDmax_right',tNDmax_right,...
    'tND_left',tND_left,'tNDmax_left',tNDmax_left, 'tNDmax', tNDmax)  ;
params  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;
% params  = bads(objFun, startValues, lBounds, uBounds, [2 .01 0 -.2], [70 1 .5 .2]) ;

[NLL, f] = GetExtremaNLL(params, d,0, 3, ...
    'k',k,'tND_right',tND_right,'tNDmax_right',tNDmax_right,...
    'tND_left',tND_left,'tNDmax_left',tNDmax_left, 'tNDmax', tNDmax) ;

p4p = params ;
p4p(4) = 0 ;

[~,f4p] = GetExtremaNLL(p4p, d,0, 1, ...
    'k',k,'tND_right',tND_right,'tNDmax_right',tNDmax_right,...
    'tND_left',tND_left,'tNDmax_left',tNDmax_left, 'tNDmax', tNDmax) ;

f.c4p_pCorrect = f.c4p(f.c4p>=0) ;
f.predChoice_pCorrect = mean([fliplr(1-f.predChoice(1:end/2)); f.predChoice( (end/2)+1:end)]) ;
f.predRTmean_pCorrect = mean([fliplr(f.predRTmean(1:end/2)); f.predRTmean( (end/2)+1:end)]) ;





end