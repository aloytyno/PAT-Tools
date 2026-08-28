function [sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir)
%IMPORTNIRSPECTRA Read spectra and extract wavelength-tagged variables.
%   [SAMPLEMETADATA, WAVELENGTHS, SPECTRA] = IMPORTNIRSPECTRA(DATADIR)
%   reads SPECTRADATA.CSV from DATADIR, identifies spectral columns using
%   the WL_ naming convention, and separates metadata from the numeric
%   spectral matrix.
%
%   Input
%   -----
%   DATADIR : string or char
%       Folder containing the raw CSV files.
%
%   Outputs
%   -------
%   SAMPLEMETADATA : table
%       Non-spectral columns from the spectra file.
%   WAVELENGTHS : numeric vector
%       Wavelength values parsed from the WL_ variable names.
%   SPECTRA : double matrix
%       Numeric spectra matrix with one row per sample.
%
%   Example
%   -------
%   [sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir);

filePath = fullfile(dataDir, "spectraData.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["SampleID", "BatchID"], "string");
spectraTable = readtable(filePath, opts);

variableNames = string(spectraTable.Properties.VariableNames);
% Spectral variables are encoded as WL_<wavelength>.
isSpectralVariable = startsWith(variableNames, "WL_");
wavelengths = str2double(extractAfter(variableNames(isSpectralVariable), "WL_"))';
spectra = spectraTable{:, isSpectralVariable};
sampleMetadata = spectraTable(:, ~isSpectralVariable);
end
