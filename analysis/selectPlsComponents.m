function [numComponents, rmseByComponent] = selectPlsComponents(xTrain, yTrain, maxComponents, numFolds)
%SELECTPLSCOMPONENTS Choose a PLS model size by cross-validation.
%   NUMCOMPONENTS = SELECTPLSCOMPONENTS(XTRAIN, YTRAIN) fits a sequence of
%   PLS models and returns the number of components with the lowest
%   cross-validated RMSE.
%
%   [NUMCOMPONENTS, RMSEBYCOMPONENT] = SELECTPLSCOMPONENTS(XTRAIN, YTRAIN,
%   MAXCOMPONENTS, NUMFOLDS) also returns the full RMSE profile across the
%   tested component counts.
%
%   Inputs
%   ------
%   XTRAIN : double matrix
%       Predictor matrix for calibration samples.
%   YTRAIN : double vector
%       Response values aligned with XTRAIN.
%   MAXCOMPONENTS : scalar integer, optional
%       Upper bound for the number of latent variables to test.
%   NUMFOLDS : scalar integer, optional
%       Number of cross-validation folds.
%
%   Outputs
%   -------
%   NUMCOMPONENTS : scalar integer
%       Component count with the minimum cross-validated RMSE.
%   RMSEBYCOMPONENT : vector
%       RMSE for each tested component count.
%
%   Example
%   -------
%   [nComp, rmseProfile] = selectPlsComponents(xTrain, yTrain, 8, 5);

if nargin < 3 || isempty(maxComponents)
    maxComponents = min(10, size(xTrain, 2) - 1);
end
if nargin < 4 || isempty(numFolds)
    numFolds = 5;
end

maxComponents = max(1, min(maxComponents, size(xTrain, 1) - 1));
% PLSREGRESS returns mean squared error for each model size, including the
% zero-component intercept-only model in the first column.
[~, ~, ~, ~, ~, ~, mse] = plsregress(xTrain, yTrain, maxComponents, "CV", numFolds);
rmseByComponent = sqrt(squeeze(mse(2, 2:end)));
[~, numComponents] = min(rmseByComponent);
end
