function referenceData = importReferenceData(dataDir)
filePath = fullfile(dataDir, "referenceData.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["SampleID", "BatchID", "LabStatus"], "string");
referenceData = readtable(filePath, opts);
end
