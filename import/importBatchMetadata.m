function batchMetadata = importBatchMetadata(dataDir)
filePath = fullfile(dataDir, "batchMetadata.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["BatchID", "Campaign", "Abnormality"], "string");
batchMetadata = readtable(filePath, opts);
batchMetadata.OperatorShift = categorical(batchMetadata.OperatorShift);
end
