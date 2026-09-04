classdef TestPatFunctionSuite < matlab.unittest.TestCase
    % End-to-end unit tests for the reusable PAT helper functions.

    properties
        % Project root used to build absolute paths during the test run.
        ProjectRoot
        % Per-test temporary folder for import fixtures and other scratch data.
        TempRoot
        % Paths added for the duration of each test method.
        AddedPaths
    end

    methods (TestMethodSetup)
        function setupEnvironment(testCase)
            % Keep each test self-contained by rebuilding the path and a
            % fresh temporary workspace for every method.
            testsDir = fileparts(mfilename("fullpath"));
            testCase.ProjectRoot = fileparts(testsDir);
            testCase.AddedPaths = { ...
                testCase.ProjectRoot, ...
                fullfile(testCase.ProjectRoot, "analysis"), ...
                fullfile(testCase.ProjectRoot, "preprocessing"), ...
                fullfile(testCase.ProjectRoot, "import")};
            addpath(testCase.AddedPaths{:});

            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
        end
    end

    methods (TestMethodTeardown)
        function teardownEnvironment(testCase)
            % Remove temporary path changes so tests do not leak state into
            % later methods or the interactive MATLAB session.
            if ~isempty(testCase.AddedPaths)
                rmpath(testCase.AddedPaths{:});
            end

            % Delete temporary files created by fixture helpers.
            if ~isempty(testCase.TempRoot) && isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, "s");
            end
        end
    end

    methods (Test)
        function testImportNirSpectraParsesWavelengthsAndMetadata(testCase)
            dataDir = createValidImportFixture(testCase.TempRoot);

            [sampleMetadata, wavelengths, spectra] = importNirSpectra(dataDir);

            % The fixture intentionally stores the spectral rows in a
            % non-sorted sample order so row preservation is testable.
            testCase.verifyEqual(wavelengths, [1000; 1004; 1008]);
            testCase.verifyEqual(string(sampleMetadata.SampleID), ["S2"; "S1"; "S3"]);
            testCase.verifyFalse(any(startsWith(string(sampleMetadata.Properties.VariableNames), "WL_")));
            testCase.verifySize(spectra, [3 3]);
            testCase.verifyEqual(spectra(1, :), [0.11 0.12 0.13], AbsTol = 1e-12);
        end

        
        function testImportNirSpectraThrowsForMissingFile(testCase)
            verifyThrowsAny(testCase, @() importNirSpectra(fullfile(testCase.TempRoot, "missing")));
        end

        function testImportReferenceDataReturnsTypedTable(testCase)
            dataDir = createValidImportFixture(testCase.TempRoot);

            referenceData = importReferenceData(dataDir);

            testCase.verifyClass(referenceData.SampleID, "string");
            testCase.verifyClass(referenceData.BatchID, "string");
            testCase.verifyEqual(referenceData.AssayPct(2), 100.2, AbsTol = 1e-12);
            testCase.verifyEqual(referenceData.LabStatus(1), "Released");
        end

        function testImportReferenceDataThrowsForMissingFile(testCase)
            verifyThrowsAny(testCase, @() importReferenceData(fullfile(testCase.TempRoot, "missing")));
        end

        function testImportBatchMetadataConvertsOperatorShift(testCase)
            dataDir = createValidImportFixture(testCase.TempRoot);

            batchMetadata = importBatchMetadata(dataDir);

            testCase.verifyClass(batchMetadata.BatchID, "string");
            testCase.verifyClass(batchMetadata.OperatorShift, "categorical");
            testCase.verifyEqual(string(batchMetadata.Abnormality), ["Normal"; "ScatterShift"]);
        end

        function testImportBatchMetadataThrowsForMissingFile(testCase)
            verifyThrowsAny(testCase, @() importBatchMetadata(fullfile(testCase.TempRoot, "missing")));
        end

        function testBuildPatDatasetJoinsDataAndRestoresRowOrder(testCase)
            dataDir = createValidImportFixture(testCase.TempRoot);

            [sampleData, wavelengths, spectra] = buildPatDataset(dataDir);

            % The joined table should follow the original spectra row order,
            % not the reference-data order used in the second fixture file.
            testCase.verifyEqual(string(sampleData.SampleID), ["S2"; "S1"; "S3"]);
            testCase.verifyEqual(string(sampleData.BatchID), ["B2"; "B1"; "B2"]);
            testCase.verifyEqual(sampleData.AssayPct, [100.2; 99.1; 98.9], AbsTol = 1e-12);
            testCase.verifyEqual(wavelengths, [1000; 1004; 1008]);
            testCase.verifyEqual(spectra(:, 1), [0.11; 0.21; 0.31], AbsTol = 1e-12);
        end

        function testBuildPatDatasetRejectsDuplicateSampleIds(testCase)
            dataDir = createDuplicateSampleFixture(testCase.TempRoot);

            testCase.verifyError(@() buildPatDataset(dataDir), "PATDemo:DuplicateSamples");
        end

        function testApplySnvNormalizesRowsAndHandlesFlatSpectrum(testCase)
            x = [1 2 3; 10 10 10];

            xSnv = applySnv(x);

            testCase.verifyEqual(mean(xSnv(1, :)), 0, AbsTol = 1e-12);
            testCase.verifyEqual(std(xSnv(1, :), 0, 2), 1, AbsTol = 1e-12);
            testCase.verifyEqual(xSnv(2, :), [0 0 0], AbsTol = 1e-12);
        end

        function testApplySavitzkyGolayPreservesSizeAndRejectsBadInputs(testCase)
            wavelengths = (1000:4:1036)';
            x = repmat(5, 2, numel(wavelengths));

            xProcessed = applySavitzkyGolay(x, 2, 5, 1, wavelengths);

            % A flat spectrum should stay numerically near zero after the
            % derivative step, while matrix shape stays unchanged.
            testCase.verifySize(xProcessed, size(x));
            testCase.verifyLessThan(max(abs(xProcessed), [], "all"), 1e-10);
            testCase.verifyError(@() applySavitzkyGolay(x, 2, 4, 1, wavelengths), "PATDemo:EvenFrameLength");
            testCase.verifyError(@() applySavitzkyGolay(x, 2, 5, 1, wavelengths(1:end-1)), "PATDemo:WavelengthMismatch");
        end

        function testMeanCenterDataUsesTrainingCenter(testCase)
            xTrain = [1 2; 3 4];
            xApply = [5 7];

            [xTrainCentered, xApplyCentered, center] = meanCenterData(xTrain, xApply);

            testCase.verifyEqual(center, [2 3], AbsTol = 1e-12);
            testCase.verifyEqual(xTrainCentered, [-1 -1; 1 1], AbsTol = 1e-12);
            testCase.verifyEqual(xApplyCentered, [3 4], AbsTol = 1e-12);

            [xTrainCenteredDefault, xApplyCenteredDefault] = meanCenterData(xTrain);
            testCase.verifyEqual(xTrainCenteredDefault, xApplyCenteredDefault, AbsTol = 1e-12);
        end

        function testRunPreprocessingPipelineBuildsReusableModel(testCase)
            wavelengths = (1000:4:1056)';
            xTrain = [sin(wavelengths / 50)'; cos(wavelengths / 55)'; sin(wavelengths / 60)' + 0.1];
            xApply = xTrain + 0.02;

            [xTrainProcessed, model] = runPreprocessingPipeline(xTrain, wavelengths);
            xApplyProcessed = runPreprocessingPipeline(xApply, wavelengths, struct(), model);

            manualTrain = applySnv(xTrain);
            manualTrain = applySavitzkyGolay(manualTrain, 2, 15, 1, wavelengths);
            manualApply = applySnv(xApply);
            manualApply = applySavitzkyGolay(manualApply, 2, 15, 1, wavelengths);

            % Compare against a manually reconstructed pipeline to ensure
            % the orchestration function is not hiding a different transform.
            testCase.verifySize(xTrainProcessed, size(xTrain));
            testCase.verifyEqual(mean(xTrainProcessed, 1), zeros(1, size(xTrain, 2)), AbsTol = 1e-10);
            testCase.verifyEqual(model.wavelengths, wavelengths);
            testCase.verifyEqual(xTrainProcessed, manualTrain - model.center, AbsTol = 1e-10);
            testCase.verifyEqual(xApplyProcessed, manualApply - model.center, AbsTol = 1e-10);
            testCase.verifyError(@() runPreprocessingPipeline(xTrain, wavelengths(1:end-1)), "PATDemo:PreprocessingSize");
        end

        function testSelectPlsComponentsReturnsBoundedProfile(testCase)
            rng(7);
            xTrain = randn(30, 4);
            yTrain = 2.5 * xTrain(:, 1) - 1.3 * xTrain(:, 2) + 0.05 * randn(30, 1);

            [numComponents, rmseByComponent] = selectPlsComponents(xTrain, yTrain, 10, 4);

            testCase.verifyGreaterThanOrEqual(numComponents, 1);
            testCase.verifyLessThanOrEqual(numComponents, size(xTrain, 2));
            testCase.verifyNumElements(rmseByComponent, size(xTrain, 2));
            testCase.verifyTrue(all(isfinite(rmseByComponent)));
        end

        function testComputePcaMonitoringStatsReturnsConsistentModel(testCase)
            rng(11);
            base = randn(40, 2);
            xTrain = [base, base(:, 1) + 0.05 * randn(40, 1), base(:, 2) - 0.05 * randn(40, 1)];

            model = computePcaMonitoringStats(xTrain, 4, 0.99);

            % Using four components on four variables should leave no
            % residual variance, so the Q limit should collapse to zero.
            testCase.verifySize(model.coeff, [4 4]);
            testCase.verifySize(model.score, [40 4]);
            testCase.verifySize(model.mu, [1 4]);
            testCase.verifyEqual(model.numComponents, 4);
            testCase.verifyGreaterThan(model.t2Limit, 0);
            testCase.verifyEqual(model.qLimit, 0, AbsTol = 1e-12);
        end

        function testApplyPcaMonitoringModelFlagsObviousOutlier(testCase)
            rng(19);
            xTrain = randn(60, 5);
            model = computePcaMonitoringStats(xTrain, 2, 0.99);
            xNew = [xTrain(1:5, :); 12 * ones(1, 5)];

            monitoring = applyPcaMonitoringModel(model, xNew);

            % Append one extreme point so the test checks that the model
            % actually raises an outlier flag rather than only returning
            % correctly shaped outputs.
            testCase.verifySize(monitoring.score, [6 2]);
            testCase.verifySize(monitoring.t2, [6 1]);
            testCase.verifySize(monitoring.qResidual, [6 1]);
            testCase.verifyTrue(monitoring.isOutlier(end));
        end
    end
end

function dataDir = createValidImportFixture(rootDir)
% Build a small but intentionally nontrivial import fixture.
% Spectra rows are ordered differently from reference rows so row-alignment
% logic in BUILDPATDATASET is exercised.
dataDir = fullfile(rootDir, "raw");
mkdir(dataDir);

spectraData = table( ...
    ["S2"; "S1"; "S3"], ...
    ["B2"; "B1"; "B2"], ...
    [10; 5; 15], ...
    [0.11; 0.21; 0.31], ...
    [0.12; 0.22; 0.32], ...
    [0.13; 0.23; 0.33], ...
    'VariableNames', ["SampleID", "BatchID", "BlendMinute", "WL_1000", "WL_1004", "WL_1008"]);

referenceData = table( ...
    ["S1"; "S2"; "S3"], ...
    ["B1"; "B2"; "B2"], ...
    [5; 10; 15], ...
    [99.1; 100.2; 98.9], ...
    [2.1; 2.3; 2.2], ...
    ["Released"; "Released"; "Review"], ...
    'VariableNames', ["SampleID", "BatchID", "BlendMinute", "AssayPct", "MoisturePct", "LabStatus"]);

batchMetadata = table( ...
    ["B1"; "B2"], ...
    ["Campaign_A"; "Campaign_B"], ...
    [220; 215], ...
    [405; 398], ...
    ["Day"; "Night"], ...
    ["Normal"; "ScatterShift"], ...
    'VariableNames', ["BatchID", "Campaign", "BlendSpeedRPM", "FeedRateKgHr", "OperatorShift", "Abnormality"]);

writetable(spectraData, fullfile(dataDir, "spectraData.csv"));
writetable(referenceData, fullfile(dataDir, "referenceData.csv"));
writetable(batchMetadata, fullfile(dataDir, "batchMetadata.csv"));
end

function dataDir = createDuplicateSampleFixture(rootDir)
% Build an import fixture that should fail duplicate-sample validation.
dataDir = fullfile(rootDir, "duplicate_raw");
mkdir(dataDir);

spectraData = table( ...
    ["S1"; "S1"], ...
    ["B1"; "B1"], ...
    [5; 10], ...
    [0.11; 0.21], ...
    [0.12; 0.22], ...
    [0.13; 0.23], ...
    'VariableNames', ["SampleID", "BatchID", "BlendMinute", "WL_1000", "WL_1004", "WL_1008"]);

referenceData = table( ...
    "S1", ...
    "B1", ...
    5, ...
    99.1, ...
    2.1, ...
    "Released", ...
    'VariableNames', ["SampleID", "BatchID", "BlendMinute", "AssayPct", "MoisturePct", "LabStatus"]);

batchMetadata = table( ...
    "B1", ...
    "Campaign_A", ...
    220, ...
    405, ...
    "Day", ...
    "Normal", ...
    'VariableNames', ["BatchID", "Campaign", "BlendSpeedRPM", "FeedRateKgHr", "OperatorShift", "Abnormality"]);

writetable(spectraData, fullfile(dataDir, "spectraData.csv"));
writetable(referenceData, fullfile(dataDir, "referenceData.csv"));
writetable(batchMetadata, fullfile(dataDir, "batchMetadata.csv"));
end

function verifyThrowsAny(testCase, fcn)
% Small helper for cases where any thrown exception is acceptable.
didThrow = false;
try
    fcn();
catch
    didThrow = true;
end
testCase.verifyTrue(didThrow);
end
