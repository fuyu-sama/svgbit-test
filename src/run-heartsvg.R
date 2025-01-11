#!/usr/bin/env Rscript

# %%
library <- function(...) suppressMessages(base::library(...))
library('dplyr')
library('data.table')
library('HEARTSVG')

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "svgbit-test")
source(fs::path(WORKDIR, "src", "utils", "read_func.R"))
source(fs::path(WORKDIR, "src", "utils", "run_func.R"))

write.dir <- fs::path(WORKDIR, "results", "HEARTSVG")
if (!fs::dir_exists(write.dir)) fs::dir_create(write.dir)

# %% DLPFC
samples <- c(
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
)
for (sample in samples) {
    reads <- read_dlpfc(sample)
    results <- run_heartsvg(reads$count_df, reads$coor_df)
    write.csv(results, fs::path(write.dir, paste0("DLPFC-", sample, ".csv")))
}

# %% iimpact
samples <- c(
    "human_breast_cancer", "human_ovarian_cancer", "human_prostate_cancer"
)
for (sample in samples) {
    reads <- read_iimpact(sample)
    results <- run_heartsvg(reads$count_df, reads$coor_df)
    write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))
}

# %% mba
for (sample in names(mba_files)) {
    result <- try({
        reads <- read_mba(sample)
        results <- run_heartsvg(reads$count_df, reads$coor_df)
        write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))
    }, silent = TRUE)
}

# %% gbm
samples <- c("21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2")
for (sample in samples) {
    reads <- read_gbm(sample)
    results <- run_heartsvg(reads$count_df, reads$coor_df)
    write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))
}

# %% stereo
sample <- "E9.5_E1S1.MOSTA"
reads <- read_stereo(sample)
results <- run_heartsvg(reads$count_df, reads$coor_df)
write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))

# %% stereo
sample <- "Mouse_brain"
reads <- read_stereo(sample)
results <- run_heartsvg(reads$count_df, reads$coor_df)
write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))

# %%
samples <- c(
    "E135A", "E135B", "E155A", "E155B", "E165A", "E165B", "E175A1", "E175A2",
    "E175B", "P0A1", "P0A2", "P0B"
)
for (sample in samples) {
    reads <- read_self(sample)
    results <- run_heartsvg(reads$count_df, reads$coor_df)
    write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))
}
