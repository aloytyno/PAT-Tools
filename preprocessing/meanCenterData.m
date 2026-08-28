function [xTrainCentered, xApplyCentered, center] = meanCenterData(xTrain, xApply)
if nargin < 2 || isempty(xApply)
    xApply = xTrain;
end

center = mean(xTrain, 1);
xTrainCentered = xTrain - center;
xApplyCentered = xApply - center;
end
