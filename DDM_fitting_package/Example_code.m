% Examples of how GMS fits diffusion and extrema models. See "help" for the
% fitting functions for more info.


addpath(genpath('/Users/gabe/Science/DTB_Fitting/DDM_fitting_package'))
load('Example_data')
%% Fit RT means and predict choice

% Get data into correct format
dStructRT = choiceMat2Struct(dataMatRT) ;

useChoiceFlag = 0 ; % Do not use choices for fitting

% Fit with a diffusion model
[paramsD, NLLd, fD] = FitAccumulationModel(dStructRT, [21 .9 .37 .37 0 7], useChoiceFlag)

% Fit with an extrema model
[paramsE, NLLe, fE] = FitExtremaModel(dStructRT, [80 .07 .4 .37 0 0], useChoiceFlag)

% Plot fits
plotBehavior_paper(dStructRT,fD,1,'m')
plotBehavior_paper(dStructRT,fE,1,'c')

%% Fit full RT distributions with collapsing bounds model (This can take a few hours)

% Diffusion
[paramsDTB, NLLdtb, fDTB] = FitDTB_CB(dStructRT, [19 1 1.2 .4 .35 .04 0 0], 'Logistic')

% Extrema
[paramsEXT, NLLext, fEXT] = FitExt_CB2(dStructRT, [100 .075 .2 1 .4 .05 0 0], 'Logistic')

%% Fit variable duration data

% Note that sCoh is in the form .512 (not 512), and stimDur is in ms.

% Get quantiles for stimulus duration and round to nearest quantile (Makes
% fitting much faster, with very little difference in fit. Also useful for plotting)
quants = quantile(dataMatVSD(:,3), 10) ; % Here I've chosen 10 quantiles.
dataMatVSD(:,3) = arbRound(quants, dataMatVSD(:,3)) ;

% Get data into the right format for fitting
dStructVSD = choiceMat2Struct_vd(dataMatVSD) ;

% Fit the data with diffusion
[paramsVSDd, NLLVSDd, fVSDd] = FitAccumulationVD(dStructVSD, [18 .8 0 0 0 0],'None') % None refers to flat bounds
% [paramsVSDd, NLLVSDd, fVSDd] = FitAccumulationVD(dStructVSD, [18 .8 0 0 0 0],'None', 'k', paramsDTB(1)) % Fit using the kappa estimated from RT fits

% Fit the data with extrema
[paramsVSDe, NLLVSDe, fVSDe] = FitExtremaVD(dStructVSD,[80 .068 0 0 0 .5], 2,'None')
% [paramsVSDe, NLLVSDe, fVSDe] = FitExtremaVD(dStructVSD,[80 .068 0 0 0 .5], 2,'None', 'k', paramsEXT(1)) 

%% Plot sensitivity as a function of stimulus duration and the model fits
figure
plotVSDdata(fVSDd, 'm')
plotVSDdata(fVSDe, 'c')


