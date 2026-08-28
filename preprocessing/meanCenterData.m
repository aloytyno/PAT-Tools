function [xTrainCentered, xApplyCentered, center] = meanCenterData(xTrain, xApply)
%MEANCENTERDATA Mean-center training data and reuse the same center.
%   [XTRAINCENTERED, XAPPLYCENTERED, CENTER] = MEANCENTERDATA(XTRAIN)
%   computes the column means of XTRAIN and subtracts them from XTRAIN.
%
%   [XTRAINCENTERED, XAPPLYCENTERED, CENTER] = MEANCENTERDATA(XTRAIN,
%   XAPPLY) applies the same center to XAPPLY so calibration and new data
%   are transformed consistently.
%
%   Inputs
%   ------
%   XTRAIN : double matrix
%       Matrix used to compute the centering vector.
%   XAPPLY : double matrix, optional
%       Matrix to center with the same training-set means.
%
%   Outputs
%   -------
%   XTRAINCENTERED : double matrix
%       Centered version of XTRAIN.
%   XAPPLYCENTERED : double matrix
%       Centered version of XAPPLY.
%   CENTER : row vector
%       Column means computed from XTRAIN.
%
%   Example
%   -------
%   [xTrainC, xValC, center] = meanCenterData(xTrain, xValidation);

if nargin < 2 || isempty(xApply)
    xApply = xTrain;
end

% Always derive the center from the calibration set.
center = mean(xTrain, 1);
xTrainCentered = xTrain - center;
xApplyCentered = xApply - center;
end
