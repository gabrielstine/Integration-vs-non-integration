function [params, NLL, f] = FitAccumulationVD(data, startValues, boundShape, varargin)

% Function to fit drift diffusion model to variable duration data via Fokker-Planck Chang-Cooper method.
%
% Created by gms. Fokker-Planck code written by az and modified by gms. Last updated
% 5/16/19.
%
% Code runs much faster if stimulus durations are rounded to nearest
% quantile, although you'll get modestly different fits.
%
% INPUTS
% 
%   %%%%%%% data (Output of choiceMat2Struct):
% 
%   1 x nCoherences struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     stimLength  (nTrials x 1, stimulus duration (in ms) for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
%
%
%   %%%%%%% startValues:
% 
%    a 1x6 array of starting parameters (if boundShape is empty or 'None', 'a' and 'd' don't do anything).
% 
%   Parameter order is:
% 
%   k         = params(1) ; (Kappa)
%   b0        = params(2) ; (Initial Bound height)
%   a         = params(3) ; (Bound collapse parameter)
%   d         = params(4) ; (Bound collapse parameter)
%   cBias     = params(5) ; (Coherence offset)
%   dBias     = params(6) ; (drift starting-point offset)
%
%
%   %%%%%%% boundShape (default is 'None'): 
% 
%     a string, specifying the shape of the bounds. Here are the possible shapes:
%     
%     'None' (Flat bounds)
%     'Linear'
%     'Quadratic'
%     'Exponential'
%     'Logistic'
%     'Hyperbolic'
%     'Step'
% 
% 
%  
%   %%%%%%% Additional inputs:
% 
%     'k', followed by a scalar: fixes kappa to the scalar value.
%
%
%
% OUTPUTS
% 
% %%%%%%% params:
% 
%     Fitted parameters, using the same order as above.
%     
%     
% %%%%%%% NLL:
% 
%     Negative-log-likelihood
%     
% 
% %%%%%%% f:
% 
%     A struct containing info for the fits. Here's an example with 8 coherences and 20 stimulus durations:
%           
%         params: [1x6 double]
%            NLL: [1x1 double]
%           sCoh: [8x1 double]         
%       stimDurs: [1×20 double]
%        nTrials: [8x20 double]       (# of trials for each stim condition)
%     predChoice: [8×20 double]       (Predicted probability rightward choice for each coherence and stim duration)
%     choiceData: [8×20 double]       (Actual probability rightward choice for each coherence and stim duration) 
%             bU: [1772×1 double]     (Upper bound)
%             bL: [1772×1 double]     (Lower bound)
%       predSens: [1×20 double]       (Predicted sensitivity (beta parameter from logistic regression) for each stimulus duration)
%       dataSens: [1x20 double]       (Actual sensitivity)
%     dataSensSE: [1x20 double]       (Standard error for dataSens)  
%             dt: .0005               (Size of time step)  
%              P: [1×1 struct]        (Results of Fokker-Planck)             
%
%
%
%                 f.P.up contains most of the useful information:
% 
%                         pdf_t: [8×1772 double]  (Probability mass above upper-bound at each time point)
%                 pdf_aboveZero: [8×1772 double]  (Probability mass >0 but not > upper-bound at each time point)
%                             p: [8×1 double]     (Total probability mass above upper-bound (up.p + lo.p should = 1 for all coherences!))
%                        mean_t: [8×1 double]     (Expectation of pdf_t (aka predicted decision times))
%                         cdf_t: [8×1772 double]  (Cumulative probability mass above upper-bound at each time point)
%                          p4VD: [8×1772 double]  (All probability mass >0 at each time point (regardless of whether its been absorbed))
%
%                 f.P.lo contains the same info but pertaining to the lower
%                 bound.
                

k = 0 ;

for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        k=varargin{i+1};
    end
        
end


options =   optimoptions(@fmincon, 'Display', 'iter', 'Maxiter', 1000, 'MaxFuneval', 2000, 'Algorithm', 'interior-point', 'UseParallel', true ) ;

lBounds(1) = 0;
lBounds(2) = 0 ;
lBounds(3) = -50 ;
lBounds(4) = -50 ;
lBounds(5) =  -.5 ;
lBounds(6) = -1 ;


uBounds(1) = 50 ;
uBounds(2) = 5 ;
uBounds(3) = 50 ;
uBounds(4) = 50 ;
uBounds(5) = .5 ;
uBounds(6) = 1 ;


pl = [5 .3   0  0  -.1 -.1 ]; 
ph = [25 2  10  6   .1  .1 ];

if isempty(startValues)
    startValues = [10 1 1 1 0 0] ;
end

if isempty(boundShape)
    boundShape = 'None' ;
end

objFun  = @(params) GetDTB_VD_FP_NLL(params, data, boundShape, 1, 'k', k) ;
% params  = fmincon(objFun,  startValues, [],[],[],[], lBounds, uBounds,[], options ) ;
params  = bads(objFun, startValues, lBounds, uBounds, pl, ph) ;

[NLL, f] = GetDTB_VD_FP_NLL(params, data, boundShape, 0, 'k', k) ;


end
