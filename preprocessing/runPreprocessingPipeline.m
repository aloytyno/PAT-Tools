function [xProcessed, model] = runPreprocessingPipeline(x, wavelengths, options, modelIn)
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
    xWorking = applySnv(xWorking);
end

xWorking = applySavitzkyGolay(xWorking, options.SgolayOrder, options.FrameLength, options.DerivativeOrder, wavelengths);

if isempty(modelIn)
    model = struct();
    model.options = options;
    model.wavelengths = wavelengths;
    if options.UseMeanCenter
        [~, xProcessed, center] = meanCenterData(xWorking, xWorking);
        model.center = center;
    else
        xProcessed = xWorking;
        model.center = zeros(1, size(xWorking, 2));
    end
else
    model = modelIn;
    if options.UseMeanCenter
        xProcessed = xWorking - model.center;
    else
        xProcessed = xWorking;
    end
end
end

function options = defaultPreprocessingOptions()
options = struct( ...
    "UseSnv", true, ...
    "SgolayOrder", 2, ...
    "FrameLength", 15, ...
    "DerivativeOrder", 1, ...
    "UseMeanCenter", true);
end

function merged = mergeOptions(defaults, overrides)
merged = defaults;
overrideFields = fieldnames(overrides);
for fieldIdx = 1:numel(overrideFields)
    merged.(overrideFields{fieldIdx}) = overrides.(overrideFields{fieldIdx});
end
end
