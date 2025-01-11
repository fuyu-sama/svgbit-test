#! /usr/bin/env Rscript

# %%
library(dplyr)
library(BayesSpace)
library(SingleCellExperiment)
library(foreach)
library(doParallel)

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "svgbit-test")
renv::activate(WORKDIR)

registerDoParallel(cores=10)

source(fs::path(WORKDIR, "src", "utils", "read_func.R"))
source(fs::path(WORKDIR, "src", "utils", "utils.R"))

perform_bayesspace <- function(sce, rank_list, step, platform, ncs) {
    foreach(method = names(rank_list)) %dopar% {
        for (begin in seq(0, step * 4, step)) {
            seed <- 10086
            set.seed(seed)
            print(paste0(
                "Now performing ", method, ": ", begin + 1, " - ", begin + step))
            sub_genes <- rank_list[[method]][(begin + 1):(begin + step)]
            sce_sub <- sce[sub_genes, ]
            sce_sub <- spatialPreprocess(
                sce_sub,
                platform = platform,
                n.PCs = 15,
                n.HVGs = step,
                skip.PCA = FALSE,
                log.normalize = FALSE,
                assay.type = "logcounts"
            )
            sce_sub <- spatialCluster(
                sce_sub,
                q = ncs,
                platform = platform,
                init.method="mclust",
                model="t"
            )
            bayesspace_result <- as.data.frame(sce_sub@colData)
            file_name <- paste0(
                sample, "-", method, "-", begin, "_", begin + step, "-", "bayesspace.csv"
            )
            write_path <- fs::path(WORKDIR, "results", "bayesspace", file_name)
            write.csv(bayesspace_result, write_path)
        }
    } -> ...
}

read_results <- function(sample) {
    rank_list <- list()

    read_df <- read.csv(
        fs::path(WORKDIR, "results", "svgbit", sample, "AI.csv"), row.names = 1
    )
    read_df <- read_df %>% arrange(desc(AI))
    rank_list$SVGbit <- rownames(read_df)

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "SOMDE", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(qval)
        rank_list$SOMDE <- read_df$g
    }

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "SpatialDE", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(qval)
        rank_list$SpatialDE <- read_df$g
    }

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "SPARK", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(adjusted_pvalue)
        rank_list$SPARK <- rownames(read_df)
    }

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "DESpace", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(FDR)
        rank_list$DESpace <- rownames(read_df)
    }

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "HEARTSVG", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(rank)
        rank_list$HEARTSVG <- read_df$gene
    }

    read_df <- try({
        read.csv(
            fs::path(WORKDIR, "results", "MERINGUE", paste0(sample, ".csv")),
            row.names = 1
        )
    })
    if (!inherits(read_df, "try-error")) {
        read_df <- read_df %>% arrange(p.adj)
        rank_list$MERINGUE <- rownames(read_df)
    }

    return(rank_list)
}

# %%
samples <- c(
    151507, 151508, 151509, 151510, 151670, 151671, 151673,
    151674, 151675, 151676
    #151669, 151672,
)
platform <- "Visium"
ncs <- 7
for (sample in samples) {
    read_files <- read_dlpfc(sample)
    sample <- paste0("DLPFC-", sample)
    rank_list <- read_results(sample)
    sce <- SingleCellExperiment(
        assays = list(
            counts = as(as.matrix(read_files$count_df), "dgCMatrix"),
            logcounts = as(as.matrix(read_files$logcount_df), "dgCMatrix")),
        colData = read_files$array_df
    )
    sce <- sce[, colSums(counts(sce)) > 0]
    perform_bayesspace(
        sce,
        rank_list,
        500,
        platform,
        ncs
    )
}

# %%
sample <- "Mouse_brain"
platform <- "ST"
ncs <- 20
read_files <- read_stereo(sample)
rank_list <- read_results(sample)
sce <- SingleCellExperiment(
    assays = list(
        counts = as(as.matrix(read_files$count_df), "dgCMatrix"),
        logcounts = as(as.matrix(read_files$logcount_df), "dgCMatrix")),
    colData = read_files$array_df
)
sce <- sce[, colSums(counts(sce)) > 0]
perform_bayesspace(
    sce,
    rank_list,
    500,
    platform,
    ncs
)
