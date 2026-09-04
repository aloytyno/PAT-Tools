function model = computePcaMonitoringStats(xTrain, numComponents, alpha)
%COMPUTEPCAMONITORINGSTATS Fit a PCA model and derive monitoring limits.
%   MODEL = COMPUTEPCAMONITORINGSTATS(XTRAIN) fits a PCA model to XTRAIN
%   using three principal components and a 95% confidence level.
%
%   MODEL = COMPUTEPCAMONITORINGSTATS(XTRAIN, NUMCOMPONENTS, ALPHA)
%   overrides the number of retained components and the confidence level
%   used for the T^2 and Q residual limits.
%
%   Inputs
%   ------
%   XTRAIN : double matrix
%       Calibration or training data with observations in rows.
%   NUMCOMPONENTS : scalar integer, optional
%       Number of principal components to retain.
%   ALPHA : scalar double, optional
%       Confidence level used to compute the monitoring limits.
%
%   Output
%   ------
%   MODEL : struct
%       Struct containing PCA coefficients, scores, latent values,
%       explained variance, centering vector, per-sample monitoring
%       statistics, and control limits.
%
%   Example
%   -------
%   pcaModel = computePcaMonitoringStats(xCalibration, 2, 0.99);

if nargin < 2 || isempty(numComponents)
    numComponents = 3;
end
if nargin < 3 || isempty(alpha)
    alpha = 0.95;
end

% Compute PCA, scores, and percent variance explained
[coeff, score, latent] = pca(xTrain);
mu = mean(xTrain, 1);
explained = 100 * latent / sum(latent);
coeff = coeff(:, 1:numComponents);
score = score(:, 1:numComponents);
% Reconstruct the training data from the retained PCA subspace.
% Compute Hotelling's T^2 for retained components
t2 = sum((score(:, 1:numComponents) .^ 2) ./ latent(1:numComponents)', 2);
xHat = score(:, 1:numComponents) * coeff(:, 1:numComponents)' + mu;
residuals = xTrain - xHat;
qResidual = sum(residuals .^ 2, 2);

numSamples = size(xTrain, 1);
% Compute T^2 control limit via F-distribution approximation
t2Limit = numComponents * (numSamples - 1) * (numSamples + 1) / ...
    (numSamples * (numSamples - numComponents)) * finv(alpha, numComponents, numSamples - numComponents);
% Compute Q residual control limit from discarded eigenvalues
qLimit = localQResidualLimit(latent, numComponents, alpha);

% Pack PCA model and monitoring statistics into struct
model = struct();
model.coeff = coeff;
model.score = score(:, 1:numComponents);
model.latent = latent;
model.mu = mu;
model.explained = explained;
model.numComponents = numComponents;
model.alpha = alpha;
model.t2 = t2;
model.qResidual = qResidual;
model.t2Limit = t2Limit;
model.qLimit = qLimit;
end

function qLimit = localQResidualLimit(latent, numComponents, alpha)
% Estimate a Q residual control limit from the discarded variance.
%LOCALQRESIDUALLIMIT - Compute Q-limit from discarded PCA eigenvalues
%
% Input arguments:
% latent        - vector of PCA eigenvalues (descending)
% numComponents - number of retained components
% alpha         - confidence level for limit (e.g., 0.95)
%
% Output arguments:
% qLimit        - scalar Q residual control limit
remainingLatent = latent(numComponents + 1:end);
% If no discarded variance, residual limit is zero
if isempty(remainingLatent)
    qLimit = 0;
    return
end

% Compute moments of the discarded eigenvalues (used in approximation)
theta1 = sum(remainingLatent);
theta2 = sum(remainingLatent .^ 2);
theta3 = sum(remainingLatent .^ 3);

% If second moment zero, variance-free residuals imply zero limit
if theta2 == 0
    qLimit = 0;
    return
end

% Calculate skewness-correction factor h0 and guard lower bound
h0 = 1 - (2 * theta1 * theta3) / (3 * theta2 ^ 2);
h0 = max(h0, 0.001);
ca = norminv(alpha);
% Assemble approximation term and raise to power 1/h0
term = ca * sqrt(2 * theta2 * h0 ^ 2) / theta1 + 1 + theta2 * h0 * (h0 - 1) / theta1 ^ 2;
qLimit = theta1 * term ^ (1 / h0);
end
