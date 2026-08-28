function [xProcessed, model] = runPreprocessingPipeline(x, wavelengths, options, modelIn)
%RUNPREPROCESSINGPIPELINE Run the configured spectral preprocessing steps.
%   [XPROCESSED, MODEL] = RUNPREPROCESSINGPIPELINE(X, WAVELENGTHS) applies
%   the default preprocessing sequence to X and returns a model struct that
%   stores reusable settings such as the mean-centering vector.
%
%   [XPROCESSED, MODEL] = RUNPREPROCESSINGPIPELINE(X, WAVELENGTHS, OPTIONS)
%   overrides the default preprocessing settings.
%
%   [XPROCESSED, MODEL] = RUNPREPROCESSINGPIPELINE(X, WAVELENGTHS, OPTIONS,
%   MODELIN) applies the pipeline using an existing preprocessing model,
%   which is typically used when transforming validation or future data
%   with the same calibration-set center.
%
%   Inputs
%   ------
%   X : double matrix
%       Spectral matrix with one sample per row.
%   WAVELENGTHS : column vector
%       Wavelengths aligned with the columns of X.
%   OPTIONS : struct, optional
%       Overrides for UseSnv, SgolayOrder, FrameLength, DerivativeOrder,
%       and UseMeanCenter.
%   MODELIN : struct, optional
%       Existing preprocessing model returned by an earlier calibration
%       call to this function.
%
%   Outputs
%   -------
%   XPROCESSED : double matrix
%       Preprocessed spectral matrix.
%   MODEL : struct
%       Struct containing the resolved options, wavelength vector, and
%       centering information for reuse.
%
%   Example
%   -------
%   [xTrainPrep, prepModel] = runPreprocessingPipeline(xTrain, wavelengths);
%   xValPrep = runPreprocessingPipeline(xValidation, wavelengths, ...
%       struct(), prepModel);

if nargin < 3 || isempty(options)
    options = defaultPreprocessingOptions();
end
if nargin < 4
    modelIn = [];
end

if size(x, 2) ~= numel(wavelengths)
    error("PATDemo:PreprocessingSize", "Input spectra must align with the wavelength vector.");
end

options = mergeOptions(defaultPreprocessingOptions(), options);
xWorking = x;

if options.UseSnv
    % SNV compensates for row-wise scatter and offset differences.
    xWorking = applySnv(xWorking);
end

xWorking = applySavitzkyGolay(xWorking, options.SgolayOrder, options.FrameLength, options.DerivativeOrder, wavelengths);

if isempty(modelIn)
    model = struct();
    model.options = options;
    model.wavelengths = wavelengths;
    if options.UseMeanCenter
        % Store the calibration-set center so later data can be aligned to
        % the same reference point.
        [~, xProcessed, center] = meanCenterData(xWorking, xWorking);
        model.center = center;
    else
        xProcessed = xWorking;
        model.center = zeros(1, size(xWorking, 2));
    end
else
    model = modelIn;
    if options.UseMeanCenter
        % Reapply the original calibration-set center to new data.
        xProcessed = xWorking - model.center;
    else
        xProcessed = xWorking;
    end
end
end

function options = defaultPreprocessingOptions()
% Default settings chosen for the NIR PAT demo workflow.
options = struct( ...
    "UseSnv", true, ...
    "SgolayOrder", 2, ...
    "FrameLength", 15, ...
    "DerivativeOrder", 1, ...
    "UseMeanCenter", true);
end

function merged = mergeOptions(defaults, overrides)
% Merge user overrides into the default preprocessing settings.
merged = defaults;
overrideFields = fieldnames(overrides);
for fieldIdx = 1:numel(overrideFields)
    merged.(overrideFields{fieldIdx}) = overrides.(overrideFields{fieldIdx});
end
end
