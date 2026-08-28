function referenceData = importReferenceData(dataDir)
%IMPORTREFERENCEDATA Read sample-level laboratory reference values.
%   REFERENCEDATA = IMPORTREFERENCEDATA(DATADIR) reads REFERENCEDATA.CSV
%   from DATADIR and returns the imported table with identifier columns
%   stored as strings.
%
%   Input
%   -----
%   DATADIR : string or char
%       Folder containing the raw CSV files.
%
%   Output
%   ------
%   REFERENCEDATA : table
%       Table of sample identifiers, batch identifiers, and reference
%       values such as assay and moisture.
%
%   Example
%   -------
%   referenceData = importReferenceData(dataDir);

filePath = fullfile(dataDir, "referenceData.csv");
opts = detectImportOptions(filePath);
opts = setvartype(opts, ["SampleID", "BatchID", "LabStatus"], "string");
referenceData = readtable(filePath, opts);
end
