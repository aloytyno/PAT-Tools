# PAT-Tools Reference

This README documents only the MATLAB functions in the `analysis`, `preprocessing`, and `import` subfolders.

## `analysis`

### `selectPlsComponents`
Selects the number of PLS components by running cross-validated `plsregress` on calibration data and choosing the component count with the lowest RMSE.

Inputs:
- `xTrain`: calibration predictor matrix
- `yTrain`: calibration response vector
- `maxComponents`: optional upper bound on the number of latent variables
- `numFolds`: optional number of cross-validation folds

Outputs:
- `numComponents`: selected component count
- `rmseByComponent`: RMSE value for each tested component count

### `computePcaMonitoringStats`
Builds a PCA monitoring model from training spectra and computes model-space statistics and control limits.

Inputs:
- `xTrain`: training matrix used to fit the PCA model
- `numComponents`: optional number of retained principal components
- `alpha`: optional confidence level used for T^2 and Q limits

Outputs:
- `model`: struct containing PCA loadings, scores, latent values, explained variance, centering vector, per-sample T^2 and Q values, and the corresponding monitoring limits

### `applyPcaMonitoringModel`
Applies a previously fitted PCA monitoring model to new data and computes projection scores, T^2, Q residuals, and outlier flags.

Inputs:
- `model`: PCA monitoring model returned by `computePcaMonitoringStats`
- `x`: new matrix to project into the PCA model

Outputs:
- `monitoring`: struct with fields `score`, `t2`, `qResidual`, and `isOutlier`

## `preprocessing`

### `applySnv`
Applies standard normal variate normalization row by row so each spectrum is centered and scaled by its own mean and standard deviation.

Inputs:
- `x`: spectral matrix

Outputs:
- `xSnv`: SNV-transformed spectral matrix

### `applySavitzkyGolay`
Applies Savitzky-Golay smoothing across spectral variables and optionally computes one or more numerical derivatives using the wavelength spacing.

Inputs:
- `x`: spectral matrix
- `polyOrder`: polynomial order for the Savitzky-Golay filter
- `frameLength`: odd smoothing window length
- `derivativeOrder`: number of derivative passes to apply after smoothing
- `wavelengths`: wavelength vector aligned with the columns of `x`

Outputs:
- `xProcessed`: smoothed and optionally differentiated spectral matrix

### `meanCenterData`
Computes the mean of the training matrix and subtracts it from both the training data and an optional second matrix.

Inputs:
- `xTrain`: matrix used to compute the column means
- `xApply`: optional matrix that should be centered with the same means

Outputs:
- `xTrainCentered`: centered version of `xTrain`
- `xApplyCentered`: centered version of `xApply`
- `center`: row vector of column means from `xTrain`

### `runPreprocessingPipeline`
Runs the configured preprocessing sequence for spectra. The current pipeline supports SNV, Savitzky-Golay smoothing and derivative processing, and optional mean centering with reusable calibration-state output.

Inputs:
- `x`: spectral matrix
- `wavelengths`: wavelength vector
- `options`: optional struct controlling SNV, filter settings, derivative order, and mean centering
- `modelIn`: optional preprocessing model whose centering values are reused when applying the pipeline to new data

Outputs:
- `xProcessed`: preprocessed spectral matrix
- `model`: struct containing the resolved options, wavelength vector, and centering information for reuse

## `import`

### `importNirSpectra`
Reads `spectraData.csv`, identifies wavelength columns by the `WL_` prefix, and separates sample metadata from the numeric spectra.

Inputs:
- `dataDir`: folder containing the raw CSV files

Outputs:
- `sampleMetadata`: non-spectral columns from the spectra file
- `wavelengths`: numeric wavelength vector parsed from the spectral variable names
- `spectra`: numeric spectral matrix

### `importReferenceData`
Reads `referenceData.csv` and returns the sample-level reference table with string-based identifier columns.

Inputs:
- `dataDir`: folder containing the raw CSV files

Outputs:
- `referenceData`: imported reference table

### `importBatchMetadata`
Reads `batchMetadata.csv`, converts identifier-like fields to strings, and converts `OperatorShift` to categorical.

Inputs:
- `dataDir`: folder containing the raw CSV files

Outputs:
- `batchMetadata`: imported batch metadata table

### `buildPatDataset`
Builds the combined PAT dataset by importing spectra, reference data, and batch metadata, joining them on sample and batch identifiers, restoring the original row order, and validating alignment.

Inputs:
- `dataDir`: folder containing the raw CSV files

Outputs:
- `sampleData`: joined sample-level metadata and reference table
- `wavelengths`: wavelength vector from the spectra import
- `spectra`: spectral matrix aligned with `sampleData`
