function [sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir)
filePath = fullfile(dataDir, "spectraData.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["SampleID", "BatchID"], "string");
spectraTable = readtable(filePath, opts);

variableNames = string(spectraTable.Properties.VariableNames);
isSpectralVariable = startsWith(variableNames, "WL_");
wavelengths = str2double(extractAfter(variableNames(isSpectralVariable), "WL_"))';
spectra = spectraTable{:, isSpectralVariable};
sampleMetadata = spectraTable(:, ~isSpectralVariable);
end
