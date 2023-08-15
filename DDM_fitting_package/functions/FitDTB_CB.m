function [params, NLL, f] = FitDTB_CB(data, startValues, boundShape, varargin)

% Function to fit a diffusion to bound (DTB) model with collapsing bounds to full, choice-conditioned RT distributions.
% Minimizes the function GetDTB_FP_NLL.
% 
%   [params, NLL, f] = FitDTB_CB(data, startValues, boundShape, varargin)
% 
% 
%   INPUTS:
% 
% 
% data takes the following form:
% 
%   data = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
% 
%     This data struct is the output of choiceMat2Struct
% 
% 
% 
%
%   startValues is a 1x8 vector of starting parameters.
% 
%   Parameter order is:
% 
% 
%         k       = params(1) ; Kappa
%         b0      = params(2) ; Start bound value
%         a       = params(3) ; Collapse parameter
%         d       = params(4) ; Collapse parameter
%         tNDmean = params(5) ; Nondecision time mean
%         tNDstd  = params(6) ; Nondecision time std
%         offset  = params(7) ; Bias parameter (coherence offset)
%         v       = params(8) ; Scaling parameter for variance. Fixed to be zero (ignored in the fits unless you change the code)
% 
% 
% 
% 
% 
%   boundShape is a string that specifies the function for the collapsing bound. Can take the following forms (we typically use logistic):
%     
%             case 'Linear'
%                 Bup(s)=B0-a*(t(s)-d);
%             case 'Quadratic'
%                 Bup(s)=B0-a*(t(s)-d).^2;
%             case 'Exponential'
%                 Bup(s)=B0*exp(-a*(t(s)-d));
%             case 'Logistic'
%                 Bup = B0*1./(1+exp(a*(t-d)));
%             case 'Hyperbolic'
%                 Bup(s) = B0*1./(1+(a*(t(s)-d)));
%             case 'Step'
%                 Bup(s)=B0*1e-3;
%             case {'None','Deadline'}
%                 % Just use b0
% 
%                 
%     Additional inputs:
%     
%         'tND'       -- fixes the nondecision time mean
%         'tNDmax'    -- sets a maximum value for the nondecision time mean
%         'tNDstd'    -- fixes the nondecision time std
%         'tNDstdMax' -- sets a maximum value for the nondecision time std
% 
% 
%         OUTPUTS
%         
%     params is the fitted parameter values (1x8)
%     NLL is the negative-log-likelihood
%     f is a structure that contains the fits
% 
%     Written by gms. Fokker-Planck code written by az and modified by gms.
%     Last updated 5/16/19


%%
% tic


options =   optimoptions(@fmincon, 'Display',   'iter',   'Maxiter',   1000,   'MaxFuneval',   2000,   'Algorithm',   'interior-point', 'UseParallel', true ) ;

lBounds(1) = 0 ;
lBounds(2) = 0 ;
lBounds(3) = 0 ;
lBounds(4) = 0 ;
lBounds(5) =  .05 ;
lBounds(6) =  .00001 ;
lBounds(7) = -.3 ;
lBounds(8) = -.1 ;

uBounds(1) = 40 ;
uBounds(2) = 5 ;
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



% pl = [2 .1 -10 -10 .1 .01 -.1 0]; 
% ph = [20 2  10  10 .5 .1   .1 1];

% startValues = [12.1800    1.5200    2.2200    0.5100    0.2100    0.0900   0.0100] ;
% 
% boundShape = 'Logistic' ;

objFun  = @(params) GetDTB_FP_NLL(params, data, boundShape, 0, 1) ;
% params  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;
% params  = bads(objFun, startValues, lBounds, uBounds, pl, ph) ;
params  = bads(objFun, startValues, lBounds, uBounds) ;

[NLL, f] = GetDTB_FP_NLL(params, data, boundShape, 0, 0) ;

% k       = params(1) ;
% b0      = params(2) ;
% a       = params(3) ;
% d       = params(4) ;
% tNDmean = params(5) ;
% tNDstd  = params(6) ;
% offset  = params(7) ;

% toc
