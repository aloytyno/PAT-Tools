function monitoring = applyPcaMonitoringModel(model, x)
numComponents = model.numComponents;
xCentered = x - model.mu;
scores = xCentered * model.coeff(:, 1:numComponents);
xHat = scores * model.coeff(:, 1:numComponents)' + model.mu;
residuals = x - xHat;

monitoring = struct();
monitoring.score = scores;
monitoring.t2 = sum((scores .^ 2) ./ model.latent(1:numComponents)', 2);
monitoring.qResidual = sum(residuals .^ 2, 2);
monitoring.isOutlier = monitoring.t2 > model.t2Limit | monitoring.qResidual > model.qLimit;
end
