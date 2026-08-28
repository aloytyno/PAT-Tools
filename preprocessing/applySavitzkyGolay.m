function xProcessed = applySavitzkyGolay(x, polyOrder, frameLength, derivativeOrder, wavelengths)
%APPLYSAVITZKYGOLAY Smooth spectra and optionally compute derivatives.
%   XPROCESSED = APPLYSAVITZKYGOLAY(X, POLYORDER, FRAMELENGTH,
%   DERIVATIVEORDER, WAVELENGTHS) applies Savitzky-Golay smoothing across
%   the spectral axis and then performs the requested number of numerical
%   derivative passes using the wavelength spacing.
%
%   Inputs
%   ------
%   X : double matrix
%       Spectral data with one sample per row.
%   POLYORDER : scalar integer
%       Polynomial order used in the Savitzky-Golay filter.
%   FRAMELENGTH : scalar integer
%       Odd window length used for smoothing.
%   DERIVATIVEORDER : scalar integer
%       Number of derivative passes to apply after smoothing.
%   WAVELENGTHS : column vector
%       Wavelength vector corresponding to the columns of X.
%
%   Output
%   ------
%   XPROCESSED : double matrix
%       Smoothed and optionally differentiated spectra.
%
%   Example
%   -------
%   xSg = applySavitzkyGolay(xSnv, 2, 15, 1, wavelengths);

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

% Smooth along the variable dimension so each row remains one spectrum.
xProcessed = smoothdata(x, 2, "sgolay", frameLength, "Degree", polyOrder);
spacing = mean(diff(wavelengths));
for derivativeIdx = 1:derivativeOrder
    % Apply simple finite differences after smoothing to obtain
    % derivative-like features while preserving matrix size.
    xProcessed = localDifferentiate(xProcessed, spacing);
end
end

function derivative = localDifferentiate(x, spacing)
% Use centered differences internally and one-sided differences at edges.
derivative = zeros(size(x));
derivative(:, 2:end-1) = (x(:, 3:end) - x(:, 1:end-2)) / (2 * spacing);
derivative(:, 1) = (x(:, 2) - x(:, 1)) / spacing;
derivative(:, end) = (x(:, end) - x(:, end-1)) / spacing;
end
