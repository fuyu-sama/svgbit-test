#! /usr/bin/env Rscript

# %%
WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/svgbit_local_test/")
renv::activate(WORKDIR)

source(paste0(WORKDIR, "src/utils/read_func.R"))
source(paste0(WORKDIR, "src/utils/run_func.R"))

idx_full <- c(
    "E135A", "E135B",
    "E155A", "E155B",
    "E165A", "E165B",
    "E175A1", "E175A2",
    "P0A1", "P0A2"
)

# %%
time.df <- t(data.frame(row.names = c("sample", "start", "end", "total")))
for (idx in idx_full) {
    read.list <- read_self(idx)
    start_time <- Sys.time()
    run_spark(read.list$count_df, read.list$coor_df)
    end_time <- Sys.time()
    t = as.numeric(end_time) - as.numeric(start_time)
    time.df <- rbind(
        time.df, c(idx, as.character(start_time), as.character(end_time), t))
}
write.csv(time.df, "time4.txt")

# %%
idx <- "151673"
reads <- read_dlpfc(idx)
start_time <- Sys.time()
results <- run_spark(reads$count_df, reads$coor_df)
end_time <- Sys.time()
t = as.numeric(end_time) - as.numeric(start_time)
time.df <- rbind(
    time.df, c(idx, as.character(start_time), as.character(end_time), t))

# %%
idx <- "Mouse_brain"
reads <- read_stereo(idx)
start_time <- Sys.time()
results <- run_spark(reads$count_df, reads$coor_df)
end_time <- Sys.time()
t = as.numeric(end_time) - as.numeric(start_time)
time.df <- rbind(
    time.df, c(idx, as.character(start_time), as.character(end_time), t))
