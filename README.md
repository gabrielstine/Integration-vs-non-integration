# Integration-vs-non-integration
Contains raw data and model fitting code from the paper, "Differentiating between integration and non-integration strategies in perceptual decision making".

The .mat file associated with each subject contais three variables -- dataFR, dataVSD, and dataSpeeded -- which correspond to three
trial types each subject performed. Each variable is an nTrials x 3 matrix with the following format: 

[signed_motion_coherence     choice(1==right, 0==left)   RT or stimulus duration].

Code for fitting the integration and extrema detection models, and example use cases, can be found in the DTB fitting package. Some of the code uses "BADS" for optimization, which can be found here: https://github.com/acerbilab/bads.
Feel free to email me at gabriel.stine@columbia.edu with any questions.
