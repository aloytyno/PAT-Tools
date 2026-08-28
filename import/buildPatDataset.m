function [sampleData, wavelengths, spectra] = buildPatDataset(dataDir)
[sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir);
referenceData = importReferenceData(dataDir);
batchMetadata = importBatchMetadata(dataDir);

sampleMetadata.RowOrder = (1:height(sampleMetadata))';
sampleData = innerjoin(sampleMetadata, referenceData(:, ["SampleID", "AssayPct", "MoisturePct", "LabStatus"]), ...
    "Keys", "SampleID");
sampleData = innerjoin(sampleData, batchMetadata, "Keys", "BatchID");
sampleData = sortrows(sampleData, "RowOrder");
spectra = spectra(sampleData.RowOrder, :);
sampleData = removevars(sampleData, "RowOrder");

if height(sampleData) ~= size(spectra, 1)
    error("PATDemo:ImportAlignment", "Imported sample metadata and spectra are misaligned.");
end
if numel(unique(sampleData.SampleID)) ~= height(sampleData)
    error("PATDemo:DuplicateSamples", "Sample identifiers must be unique after import.");
end
end
