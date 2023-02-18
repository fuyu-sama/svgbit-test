#! /usr/bin/env Rscript

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#  Codes are far away from bugs with Buddha's bless
#

# %%
library(dplyr)
library(BayesSpace)
library(SingleCellExperiment)
library(foreach)
library(doParallel)

WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/svgbit_local_test/")
renv::activate(WORKDIR)

registerDoParallel(cores=5)

source(paste0(WORKDIR, "src/utils/read_func.R"))
source(paste0(WORKDIR, "src/utils/utils.R"))

perform_bayesspace <- function(sce, rank_list, step, platform, ncs) {
    for (method in names(rank_list)) {
    #foreach(method = names(rank_list)) %dopar% {
        #for (begin in seq(0, step * 4, step)) {
        foreach(begin = seq(0, step * 4, step)) %dopar% {
            seed <- 42
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
            if (exists("seed")) {
                write_path <- paste0(
                    WORKDIR, "results/bayesspace/seed-", seed, "/",
                    sample, "-", method, "-", begin, "_", begin + step,
                    "-", "bayesspace.csv"
                )
            } else {
                write_path <- paste0(
                    WORKDIR, "results/bayesspace/",
                    sample, "-", method, "-", begin, "_", begin + step,
                    "-", "bayesspace.csv"
                )
            }
            write.csv(bayesspace_result, write_path)
        } -> ...
    }
}

read_results <- function(sample) {
    rank_list <- list()

    read_df <- read.csv(paste0(WORKDIR, "results/", sample, "/AI.csv"), row.names = 1)
    read_df <- read_df %>% arrange(desc(AI))
    rank_list$SVGbit <- rownames(read_df)

    try({
        read_df <- read.csv(
            paste0(WORKDIR, "results/somde/", sample, "-somde.csv"),
            row.names = 1
        )
        read_df <- read_df %>% arrange(qval)
        rank_list$SOMDE <- read_df$g
    })

    try({
        read_df <- read.csv(
            paste0(WORKDIR, "results/spatialde/", sample, "-spatialde.csv"),
            row.names = 1
        )
        read_df <- read_df %>% arrange(qval)
        rank_list$SpatialDE <- read_df$g
    })

    try({
        read_df <- read.csv(
            paste0(WORKDIR, "results/spark/", sample, "-spark.csv"),
            row.names = 1
        )
        read_df <- read_df %>% arrange(adjusted_pvalue)
        rank_list$SPARK <- rownames(read_df)
    })

    return(rank_list)
}

# %%
sample <- 151673
platform <- "Visium"
ncs <- 7
read_files <- read_dlpfc(sample)
sample <- paste0("DLPFC-", sample)
rank_list <- read_results(sample)
sce <- SingleCellExperiment(
    assays = list(
        counts = as(as.matrix(read_files$count_df), "dgCMatrix"),
        logcounts = as(as.matrix(read_files$logcounts_df), "dgCMatrix")),
    colData = read_files$array_df
)
sce <- sce[, colSums(counts(sce)) > 0]
for (step in c(500)) {
    perform_bayesspace(
        sce,
        rank_list,
        step,
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
array_df <- split_coor(colnames(read_files$count_df))
colnames(array_df) <- c("row", "col")
sce <- SingleCellExperiment(
    assays = list(
        counts = as(as.matrix(read_files$count_df), "dgCMatrix"),
        logcounts = as(as.matrix(read_files$logcounts_df), "dgCMatrix")),
    colData = array_df
)
sce <- sce[, colSums(counts(sce)) > 0]
for (step in c(500)) {
    perform_bayesspace(
        sce,
        rank_list,
        step,
        platform,
        ncs
    )
}
