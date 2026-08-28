function xProcessed = applySavitzkyGolay(x, polyOrder, frameLength, derivativeOrder, wavelengths)
arguments
    x double
    polyOrder (1,1) double {mustBeInteger, mustBeNonnegative}
    frameLength (1,1) double {mustBeInteger, mustBePositive}
    derivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative}
    wavelengths (:,1) double
end

if mod(frameLength, 2) == 0
    error("PATDemo:EvenFrameLength", "Savitzky-Golay frame length must be odd.");
end
if frameLength <= polyOrder
    error("PATDemo:InvalidFrameLength", "Frame length must be greater than polynomial order.");
end
if size(x, 2) ~= numel(wavelengths)
    error("PATDemo:WavelengthMismatch", "The number of spectral variables must match the wavelength vector.");
end

xProcessed = smoothdata(x, 2, "sgolay", frameLength, "Degree", polyOrder);
spacing = mean(diff(wavelengths));
for derivativeIdx = 1:derivativeOrder
    xProcessed = localDifferentiate(xProcessed, spacing);
end
end

function derivative = localDifferentiate(x, spacing)
derivative = zeros(size(x));
derivative(:, 2:end-1) = (x(:, 3:end) - x(:, 1:end-2)) / (2 * spacing);
derivative(:, 1) = (x(:, 2) - x(:, 1)) / spacing;
derivative(:, end) = (x(:, end) - x(:, end-1)) / spacing;
end
