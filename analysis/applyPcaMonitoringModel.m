function monitoring = applyPcaMonitoringModel(model, x)
%APPLYPCAMONITORINGMODEL Apply a PCA monitoring model to new spectra.
%   MONITORING = APPLYPCAMONITORINGMODEL(MODEL, X) projects X into the PCA
%   model returned by COMPUTEPCAMONITORINGSTATS and computes scores,
%   Hotelling's T^2 values, Q residuals, and outlier flags.
%
%   Inputs
%   ------
%   MODEL : struct
%       PCA monitoring model created by COMPUTEPCAMONITORINGSTATS.
%   X : double matrix
%       New observations arranged as rows and variables as columns.
%
%   Output
%   ------
%   MONITORING : struct
%       Struct with fields SCORE, T2, QRESIDUAL, and ISOUTLIER.
%
%   Example
%   -------
%   monitoring = applyPcaMonitoringModel(model, xValidation);

numComponents = model.numComponents;
% Reuse the calibration-set mean stored in the PCA model.
xCentered = x - model.mu;
scores = xCentered * model.coeff(:, 1:numComponents);
xHat = scores * model.coeff(:, 1:numComponents)' + model.mu;
residuals = x - xHat;

monitoring = struct();
monitoring.score = scores;
% T^2 measures leverage within the retained PCA model space.
monitoring.t2 = sum((scores .^ 2) ./ model.latent(1:numComponents)', 2);
% Q residuals capture unexplained variation outside the PCA model.
monitoring.qResidual = sum(residuals .^ 2, 2);
monitoring.isOutlier = monitoring.t2 > model.t2Limit | monitoring.qResidual > model.qLimit;
end
