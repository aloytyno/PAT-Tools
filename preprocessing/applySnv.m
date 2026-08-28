function xSnv = applySnv(x)
%APPLYSNV Apply standard normal variate normalization row by row.
%   XSNV = APPLYSNV(X) mean-centers and scales each row of X using that
%   row's own standard deviation. This is commonly used for spectra when
%   scatter effects vary from sample to sample.
%
%   Input
%   -----
%   X : double matrix
%       Spectral matrix with one observation per row.
%
%   Output
%   ------
%   XSNV : double matrix
%       SNV-normalized version of X.
%
%   Example
%   -------
%   xSnv = applySnv(xRaw);

rowMean = mean(x, 2);
rowStd = std(x, 0, 2);
% Protect against division by zero for flat spectra.
rowStd(rowStd == 0) = 1;
xSnv = (x - rowMean) ./ rowStd;
end
