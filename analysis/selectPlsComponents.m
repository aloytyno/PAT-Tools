function [numComponents, rmseByComponent] = selectPlsComponents(xTrain, yTrain, maxComponents, numFolds)
if nargin < 3 || isempty(maxComponents)
    maxComponents = min(10, size(xTrain, 2) - 1);
end
if nargin < 4 || isempty(numFolds)
    numFolds = 5;
end

maxComponents = max(1, min(maxComponents, size(xTrain, 1) - 1));
[~, ~, ~, ~, ~, ~, mse] = plsregress(xTrain, yTrain, maxComponents, "CV", numFolds);
rmseByComponent = sqrt(squeeze(mse(2, 2:end)));
[~, numComponents] = min(rmseByComponent);
end
