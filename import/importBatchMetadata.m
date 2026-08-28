function batchMetadata = importBatchMetadata(dataDir)
%IMPORTBATCHMETADATA Read batch-level metadata for the PAT dataset.
%   BATCHMETADATA = IMPORTBATCHMETADATA(DATADIR) reads BATCHMETADATA.CSV
%   from DATADIR and converts key identifier fields to MATLAB strings.
%
%   Input
%   -----
%   DATADIR : string or char
%       Folder containing the raw CSV files.
%
%   Output
%   ------
%   BATCHMETADATA : table
%       Batch metadata table with OperatorShift stored as categorical.
%
%   Example
%   -------
%   batchMetadata = importBatchMetadata(dataDir);

filePath = fullfile(dataDir, "batchMetadata.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["BatchID", "Campaign", "Abnormality"], "string");
batchMetadata = readtable(filePath, opts);
% Operator shift is better treated as a discrete label than free text.
batchMetadata.OperatorShift = categorical(batchMetadata.OperatorShift);
end
