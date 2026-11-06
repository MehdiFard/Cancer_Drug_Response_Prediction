# Pipeline: from the raw data to the figures. Each rule runs one script in scripts/.
#
#   snakemake -n            show what would run
#   snakemake --cores 12    run everything that is missing or out of date
#
# Run from the project root with R and the packages listed in README.md on the PATH
# (e.g. inside the conda environment). A step reruns when its script, the R/ helpers it
# uses or its inputs change. Logs go to results/logs/.

DATASETS     = ["CCLE", "GDSC_broad", "GDSC_sanger", "PDX_all"]
CELL_LINES   = ["CCLE", "GDSC_broad", "GDSC_sanger"]
FEATURE_SETS = ["multiomics", "panel"]
CORRECTIONS  = ["uncorrected", "corrected"]
FIGURES      = ["fig1_data_and_cv", "fig2_cross_dataset", "fig3_pdx", "fig4_top_features"]

PROC   = "data/processed"
TABLES = "results/tables"
HELPERS = ["R/feature_filter.R", "R/model_utils.R"]


rule all:
    input:
        expand("results/figures/{fig}.png", fig=FIGURES),
        f"{TABLES}/summary.csv",


rule prepare_data:
    input:
        script="scripts/01_prepare_data.R",
        raw=expand("data/raw/{ds}", ds=DATASETS),
        xeva="data/raw/xeva/Xeva_PDXE_tsv.tar.gz",
    output:
        expand(PROC + "/{ds}_clean_{fs}.rds", ds=DATASETS, fs=FEATURE_SETS),
        f"{PROC}/PDX_xeva_response.rds",
    log:
        "results/logs/01_prepare_data.log",
    shell:
        "Rscript {input.script} > {log} 2>&1"


rule batch_correction:
    input:
        "scripts/02_batch_correction.R",
        "R/batch_correction.R",
        expand(PROC + "/{ds}_clean_{fs}.rds", ds=DATASETS, fs=FEATURE_SETS),
    output:
        expand(PROC + "/{ds}_corrected_{fs}.rds", ds=DATASETS, fs=FEATURE_SETS),
        f"{TABLES}/batch_effect_pca.csv",
        f"{TABLES}/batch_effect_mixing.csv",
    log:
        "results/logs/02_batch_correction.log",
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule model_comparison:
    input:
        "scripts/03_model_comparison.R",
        HELPERS,
        expand(PROC + "/{ds}_clean_multiomics.rds", ds=CELL_LINES),
    output:
        f"{TABLES}/model_comparison_cv.csv",
        f"{TABLES}/model_comparison_summary.csv",
    log:
        "results/logs/03_model_comparison.log",
    threads: workflow.cores
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule train_models:
    input:
        "scripts/04_train_models.R",
        HELPERS,
        f"{TABLES}/model_comparison_summary.csv",
        expand(PROC + "/{ds}_clean_{fs}.rds", ds=CELL_LINES, fs=FEATURE_SETS),
        expand(PROC + "/{ds}_corrected_{fs}.rds", ds=CELL_LINES, fs=FEATURE_SETS),
    output:
        [directory(f"models/{fs}_{cr}") for fs in FEATURE_SETS for cr in CORRECTIONS],
        f"{TABLES}/cv_metrics.csv",
        f"{TABLES}/cv_predictions.csv",
    log:
        "results/logs/04_train_models.log",
    threads: workflow.cores
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule cross_dataset:
    input:
        "scripts/05_cross_dataset.R",
        HELPERS,
        f"{TABLES}/cv_metrics.csv",
        expand("models/{fs}_{cr}", fs=FEATURE_SETS, cr=CORRECTIONS),
    output:
        f"{TABLES}/cross_dataset_metrics.csv",
        f"{TABLES}/cross_dataset_summary.csv",
        f"{TABLES}/cross_dataset_by_dataset.csv",
    log:
        "results/logs/05_cross_dataset.log",
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule predict_pdx:
    input:
        "scripts/06_predict_pdx.R",
        HELPERS,
        f"{TABLES}/cv_metrics.csv",
        f"{TABLES}/cross_dataset_metrics.csv",
        f"{PROC}/PDX_xeva_response.rds",
        expand("models/{fs}_{cr}", fs=FEATURE_SETS, cr=CORRECTIONS),
    output:
        f"{TABLES}/pdx_model_choice.csv",
        f"{TABLES}/pdx_predictions.csv",
        f"{TABLES}/pdx_metrics.csv",
    log:
        "results/logs/06_predict_pdx.log",
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule top_features:
    input:
        "scripts/07_top_features.R",
        HELPERS,
        "models/multiomics_corrected",
        expand(PROC + "/{ds}_corrected_multiomics.rds", ds=CELL_LINES),
    output:
        f"{TABLES}/top_features.csv",
    log:
        "results/logs/07_top_features.log",
    shell:
        "Rscript {input[0]} > {log} 2>&1"


rule figures:
    input:
        "scripts/08_figures.R",
        f"{TABLES}/batch_effect_pca.csv",
        f"{TABLES}/batch_effect_mixing.csv",
        f"{TABLES}/model_comparison_cv.csv",
        f"{TABLES}/model_comparison_summary.csv",
        f"{TABLES}/cv_metrics.csv",
        f"{TABLES}/cross_dataset_metrics.csv",
        f"{TABLES}/pdx_metrics.csv",
        f"{TABLES}/pdx_predictions.csv",
        f"{TABLES}/top_features.csv",
    output:
        expand("results/figures/{fig}.png", fig=FIGURES),
        f"{TABLES}/summary.csv",
        f"{TABLES}/cv_r2_best.csv",
    log:
        "results/logs/08_figures.log",
    shell:
        "Rscript {input[0]} > {log} 2>&1"
