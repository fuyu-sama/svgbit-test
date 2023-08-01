#! /usr/bin/env Rscript

# %%
WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/svgbit_local_test/")
renv::activate(WORKDIR)

source(paste0(WORKDIR, "src/utils/read_func.R"))
source(paste0(WORKDIR, "src/utils/run_func.R"))

# %% als
sample <- "GSM3399149_CN68_E1"
reads <- read_als(sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
rownames(count_df) <- make.unique(reads$gene_names)
results <- run_spark(count_df, coor_df)
write.csv(results, paste0(WORKDIR, "results/spark/", sample, "-spark.csv"))

# %% mba
gsm <- "GSM4459950"
sample <- "18A"
reads <- read_mba(gsm, sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, paste0(WORKDIR, "results/spark/", sample, "-spark.csv"))

# %% srereo
sample <- "E9.5_E1S1.MOSTA"
reads <- read_stereo(sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, paste0(WORKDIR, "results/spark/", sample, "-spark.csv"))

# %% srereo
sample <- "Mouse_brain"
reads <- read_stereo(sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, paste0(WORKDIR, "results/spark/", sample, "-spark.csv"))
