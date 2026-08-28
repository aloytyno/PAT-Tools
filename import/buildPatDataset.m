function [sampleData, wavelengths, spectra] = buildPatDataset(dataDir)
%BUILDPATDATASET Assemble imported PAT tables into one aligned dataset.
%   [SAMPLEDATA, WAVELENGTHS, SPECTRA] = BUILDPATDATASET(DATADIR) imports
%   the spectra, reference data, and batch metadata found in DATADIR and
%   joins them into one sample-level table.
%
%   Inputs
%   ------
%   DATADIR : string or char
%       Folder containing the raw CSV files consumed by the import layer.
%
%   Outputs
%   -------
%   SAMPLEDATA : table
%       Joined sample metadata, reference values, and batch metadata.
%   WAVELENGTHS : numeric vector
%       Wavelength vector extracted from the spectra file.
%   SPECTRA : double matrix
%       Spectral matrix aligned row-for-row with SAMPLEDATA.
%
%   Example
%   -------
%   [sampleData, wavelengths, spectra] = buildPatDataset(dataDir);

[sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir);
referenceData = importReferenceData(dataDir);
batchMetadata = importBatchMetadata(dataDir);

% Track the original spectral row order before table joins reorder rows.
sampleMetadata.RowOrder = (1:height(sampleMetadata))';
sampleData = innerjoin(sampleMetadata, referenceData(:, ["SampleID", "AssayPct", "MoisturePct", "LabStatus"]), ...
    "Keys", "SampleID");
sampleData = innerjoin(sampleData, batchMetadata, "Keys", "BatchID");
sampleData = sortrows(sampleData, "RowOrder");
spectra = spectra(sampleData.RowOrder, :);
sampleData = removevars(sampleData, "RowOrder");

% Guard against silent misalignment between metadata rows and spectra rows.
if height(sampleData) ~= size(spectra, 1)
    error("PATDemo:ImportAlignment", "Imported sample metadata and spectra are misaligned.");
end
if numel(unique(sampleData.SampleID)) ~= height(sampleData)
    error("PATDemo:DuplicateSamples", "Sample identifiers must be unique after import.");
end
end
