function [params, NLL, f] = FitAccumulationModel( d, startValues, useChoiceFlag, varargin)

% [params NLL f] = FitAccumulationModel( d, startValues, useChoiceFlag)
% Function to fit choice/mean RT data or just mean RTs with a drift
% diffusion model. This function minimizes the objective function
% GetAccumulationNLL.
% 
% d can be generated with function choiceMat2Struct. d takes the following form:
% 
%   d = 
% 
%   1 x nCoherences struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
%
%
%
%   useChoiceFlag takes the value 0, 1, or 3 (don't ask)
%       if 0, calculate NLL using only meanRTs. If you want to fit RTs and
%       predict choice, use this one.
%       if 1, calculate NLL using both RTs and choice data
%
%
%   startValues is a 1x6 vector of starting parameters.
% 
%   Parameter order is:
% 
% k         = params(1) ; % Kappa
% A1        = params(2) ; % Bound
% tND_right = params(3) ; % non-decision time, righward
% tND_left  = params(4) ; % non-decision time, leftward
% offset    = params(5) ; % Bias
% v         = params(6) ; % Variance parameter (0 if no scaling)
% 
% By default, tND_left is ignored and tND_right is used as the tND for both
% choices.
%
%
% Other inputs:
% 
% You can force/constrain certain parameters when fitting. At the end of the input,
% put the parameter you'd like to constrain in quotes, followed by the value you 
% want it to be. For example, if I want to force 'k' to be '10', I would do:
%     
%     [params, NLL, f] = FitAccumulationModel(d,startValues,useChoiceFlag, 'k', 10) ;
%     
% You can also do this with 'tND_left' and tND_right'. You can also constrain the tND
% to be below a certain value. Use tNDmax_right or tNDmax_left for this.
% 
% 
% OUTPUTS
% 
%     params is the list of fitted parameters in the order shown above.
%     
%     NLL is the negative log-likelihood of the CHOICE data (See NLL2 for
%     the returned value of the objective function)
%     
%     f is a structure that contains more info, like the predicted psychometric 
%     and chronometric functions
% 
%     Written by gms
%     Last updated 5/16/19



k = 0 ;
B = 0 ;
tND_right = 0 ;
tND_left = 0 ;
tNDmax_right = 0 ;
tNDmax_left  = 0 ;
tNDmax = 0 ;

options =   optimset( 'Display',   'off',   'Maxiter',   1000,   'MaxFuneval',   4000,   'Algorithm',   'interior-point' ) ;

lBounds(1) = 0 ;
lBounds(2) = 0 ;
lBounds(3) = 0 ;
lBounds(4) = 0 ;
lBounds(5) = -.5 ;
lBounds(6) = 0 ;


uBounds(1) = 100 ;
uBounds(2) = 20 ;
uBounds(3) = 1 ;
uBounds(4) = 1 ;
uBounds(5) = .5 ;
uBounds(6) = 30 ;

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        k=varargin{i+1};
    end
    
    if isequal(varargin{i},'B')
        B=varargin{i+1};
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

objFun  = @(params) GetAccumulationNLL(params, d, 0, useChoiceFlag, 0,...
    'k',k, 'B', B,'tND_right',tND_right,'tNDmax_right',tNDmax_right,...
    'tND_left',tND_left,'tNDmax_left',tNDmax_left, 'tNDmax', tNDmax)  ;


[params, NLL2]  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;

[NLL, f] = GetAccumulationNLL(params, d, 0, 3, 1,...
    'k',k, 'B', B, 'tND_right',tND_right,'tNDmax_right',tNDmax_right,...
    'tND_left',tND_left,'tNDmax_left',tNDmax_left, 'tNDmax', tNDmax)  ;


f.c4p_pCorrect = f.c4p(f.c4p>=0) ;
f.predChoice_pCorrect = mean([fliplr(1-f.predChoice(1:end/2)); f.predChoice( (end/2)+1:end)]) ;
f.predRTmean_pCorrect = mean([fliplr(f.predRTmean(1:end/2)); f.predRTmean( (end/2)+1:end)]) ;
f.NLL2 = NLL2 ;


end