# Host–microbiome dynamics in the MASS cohort

Analysis code and prepared inputs for the study **“Longitudinal host–microbe profiling reveals a transferable mortality signature in severe pneumonia.”** Analyses use MASS for discovery, CMAISE for external evaluation of host-response patterns and prediction models, and EARLI for comparison of virus-associated transcriptional programs.

## Repository contents

The public repository is kept concise and contains two working directories:

| Directory | Contents |
|---|---|
| `Inputs/` | Prepared analysis inputs, including cohort metadata, expression and microbial data objects, module definitions, prediction folds, frozen models, and selected intermediate results |
| `code/` | R scripts for data preparation, statistical analyses, visualization, and validation, together with shared helper functions and an analysis runner |

Analyses cover microbial profiling, host pathways and modules, HCMV-associated responses, consensus transcriptomic subtype (CTS) dynamics, and mortality prediction. Figure-based script names identify analysis components; some supplementary-figure numbers in filenames reflect earlier manuscript versions. Use the current manuscript legends and each script's inputs and outputs to identify the corresponding analysis.

`Outputs/` is a local results directory generated when scripts run and is not included in the public repository. The manuscript data package is supplied separately and contains `source_processed_data.zip`, `Source_Data.xlsx`, and `Supplementary_Data.xlsx`. These files provide processed MASS data, figure-level numerical data, and supplementary results, respectively.

## Software environment

The recorded analysis environment uses R 4.5.0. Principal dependencies include:

- Expression and enrichment: `edgeR`, `limma`, `GSVA`, `fgsea`, `clusterProfiler`, `AnnotationDbi`, and `org.Hs.eg.db`.
- Longitudinal and survival modeling: `lme4`, `lmerTest`, `nlme`, `emmeans`, `survival`, and `JM`.
- Classification and prediction: `ConsensusTranscriptomicSubtype`, `sva`, `randomForest`, `glmnet`, and `pROC`.
- Other analyses and visualization: `xCell`, `MaAsLin3`, `vegan`, `tidyverse`, `ComplexHeatmap`, and `circlize`.

Individual scripts declare their additional dependencies. Package versions used for specific analyses are described in the accompanying Supplementary Methods. Record `sessionInfo()` when reproducing an analysis.

## Running the analyses

Run commands from the repository root so that relative paths resolve correctly. Scripts start from prepared data rather than raw sequencing reads. Create local output directories before running analyses:

```sh
mkdir -p Outputs/plot_components
```

List the stages supported by the analysis runner:

```sh
Rscript code/run_current_analysis.R --list
```

Run selected stages in dependency order, for example:

```sh
Rscript code/run_current_analysis.R \
  Figure4C_FigureS9CD_joint_models.R \
  validate_CMAISE_models.R \
  Figure4C_FigureS9CD_joint_plot.R
```

Other scripts can be run directly after their required inputs or upstream results are available. For example, external application of the frozen prediction models is followed by plotting:

```sh
Rscript code/Figure6D_analysis.R
Rscript code/Figure6D_plot.R
```

Commands write numerical results and plots to local output paths and may overwrite existing results. Reuse saved model fits and prepared results when only plotting or downstream summaries need to change.

`Rscript code/run_current_analysis.R --all` runs every stage listed by the runner, including model fitting. It covers selected host-response, CTS, and model workflows rather than every script in `code/`. Check the input requirements below before using this option.

## Input requirements and analysis conventions

Input filenames containing dates identify source snapshots. Preserve sample identifiers and supplied participant-level cross-validation folds when reproducing the analyses.

Two scripts retain paths outside the public repository: `code/Figure6B_FigureS11A_analysis.R` reads saved out-of-fold predictions from the original analysis directory, and `code/validate_module_scores.R` writes a reconstruction object to an external archive directory. Supply the required prediction file and configure these paths before running those stages. The public subset does not include those external directories.

The primary MASS endpoint is 28-day all-cause mortality. The CMAISE endpoint is recorded in-hospital death within 28 days; its joint models censor live discharge at discharge and cap follow-up at day 28. The manuscript Methods and figure legends specify analysis populations, covariates, tests, and multiple-testing families.

CTS categorical analyses use the stored assignments in `Inputs/1616_meta_model.rdata` for MASS and `Inputs/260830_SRR_CTS_classification.csv` for CMAISE. Probability models use the stored class probabilities. MASS module analyses use a common score matrix, subset for each analysis, rather than repeating normalization within subgroups.

Candidate genes, modules, and expression preprocessing were prepared before prediction resampling; internal cross-validation evaluates model fitting conditional on those choices. External prediction applies frozen MASS coefficients without refitting or recalibration. Joint models that include eventual outcome in the longitudinal component describe retrospective outcome-associated trajectories.
