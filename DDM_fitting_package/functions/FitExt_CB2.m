function [params, NLL, f] = FitExt_CB2(data, startValues, boundShape, varargin)

% Function to fit an extrema detection model t full, choice-conditioned RT distributions. Minimizes the function GetExtremaColBounds_sepRT.
% 
% 
%   See help for FitDTB_CB. This function has the same parameterization.
% 
%     It's very easy for the model to end up in a strange regime, so it's important to have reasonable starting parameters.
%     Something like the following should work, but it'd be good to verify that the fit makes sense:
%     
%     [80 .075 .2 2 .4 .05 0 0]
%     
% 
%   Written by gms
%   Last updated 5/16/19

%%
% tic
options =   optimoptions(@fmincon, 'Display',   'iter',   'Maxiter',   1000,   'MaxFuneval',   2000,   'Algorithm',   'interior-point' ) ;

lBounds(1) = 0 ;
lBounds(2) = 0 ;
lBounds(3) = 0 ;
lBounds(4) = 0 ;
lBounds(5) =  .005 ;
lBounds(6) =  .0001 ;
lBounds(7) = -.3 ;
lBounds(8) = -.1 ;

uBounds(1) = 300 ;
uBounds(2) = .2 ;
uBounds(3) = 50 ;
uBounds(4) = 50 ;
uBounds(5) = 1 ;
uBounds(6) = 1 ;
uBounds(7) = .3 ;
uBounds(8) = 5 ;


for i=1:length(varargin)
    
    if isequal(varargin{i},'tND')
        tND=varargin{i+1};
        lBounds(5) = tND ;
        uBounds(5) = tND ;
        startValues(5) = tND ;
    end
    
    if isequal(varargin{i},'tNDmax')
        tNDmax=varargin{i+1};
        uBounds(5) = tNDmax ;
        startValues(5) = tNDmax - .05 ;
    end
    
    if isequal(varargin{i},'tNDstd')
        tNDstd=varargin{i+1};
        lBounds(6) = tNDstd ;
        uBounds(6) = tNDstd ;
        startValues(6) = tNDstd ;
    end
    
    if isequal(varargin{i},'tNDstdMax')
        tNDstdMax=varargin{i+1};
        uBounds(6) = tNDstdMax ;
        startValues(6) = tNDstdMax - .005 ;
    end
    
end

% pl = [20 .05   0  0 .1 .01 -.1 0]; 
% ph = [220 .1  10  6 .5 .1   .1 1];

% startValues = [12.1800    1.5200    2.2200    0.5100    0.2100    0.0900   0.0100] ;
% 
% boundShape = 'Logistic' ;

objFun  = @(params) GetExtremaColBoundsNLL_sepRT(params, data, boundShape, 0, 1) ;
% params  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;
% params  = bads(objFun, startValues, lBounds, uBounds, pl, ph) ;
params  = bads(objFun, startValues, lBounds, uBounds) ;

[NLL, f] = GetExtremaColBoundsNLL_sepRT(params, data, boundShape, 0, 0) ;

% k       = params(1) ;
% b0      = params(2) ;
% a       = params(3) ;
% d       = params(4) ;
% tNDmean = params(5) ;
% tNDstd  = params(6) ;
% offset  = params(7) ;

% toc
