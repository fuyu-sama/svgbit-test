#! /usr/bin/env Rscript

# %%
WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "svgbit-test")

source(fs::path(WORKDIR, "src", "utils", "read_func.R"))
source(fs::path(WORKDIR, "src", "utils", "run_func.R"))

write.dir <- fs::path(WORKDIR, "results", "SPARK")
if (!fs::dir_exists(write.dir)) fs::dir_create(write.dir)

# %% dlpfc
samples <- c(
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
)
for (sample in samples) {
    reads <- read_dlpfc(sample)
    count_df <- reads$count_df
    coor_df <- reads$coor_df
    results <- run_spark(count_df, coor_df)
    write.csv(results, fs::path(write.dir, paste0("DLPFC-", sample, ".csv")))
}

# %% srereo
sample <- "E9.5_E1S1.MOSTA"
reads <- read_stereo(sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))

# %% srereo
sample <- "Mouse_brain"
reads <- read_stereo(sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))

# %% mba
gsm <- "GSM4459950"
sample <- "18A"
reads <- read_mba(gsm, sample)
count_df <- reads$count_df
coor_df <- reads$coor_df
results <- run_spark(count_df, coor_df)
write.csv(results, fs::path(write.dir, paste0(sample, ".csv")))
