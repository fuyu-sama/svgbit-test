#! /usr/bin/env Rscript

# %%
library(dplyr)
library(SpatialPCA)

WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/svgbit_local_test/")
renv::activate(WORKDIR)
set.seed(42)

source(paste0(WORKDIR, "src/utils/read_func.R"))

sub_spatialpca <- function(spca_obj, shape, ncs) {
    spca_obj <- SpatialPCA_buildKernel(
        spca_obj,
        kerneltype = "gaussian",
        bandwidthtype = "SJ",
        bandwidth.set.by.user = NULL
    )
    spca_obj <- SpatialPCA_EstimateLoading(
        spca_obj,
        fast = FALSE,
        SpatialPCnum = 20
    )
    spca_obj <- SpatialPCA_SpatialPCs(spca_obj, fast = FALSE)
    clusterlabel <- walktrap_clustering(
        clusternum = ncs,
        latent_dat = spca_obj@SpatialPCs,
        knearest = round(sqrt(dim(spca_obj@SpatialPCs)[2]))
    )
    clusterlabel_refine <- refine_cluster_10x(
        clusterlabels = clusterlabel,
        location = spca_obj@location,
        shape = shape
    )
    return(clusterlabel_refine)
}

perform_spatialpca <- function(count_df, coor_df, rank_list, step, shape, ncs) {
    for (method in names(rank_list)) {
        for (begin in seq(0, step * 4, step)) {
            sub_genes <- rank_list[[method]][(begin + 1):(begin + step)]
            count_sub <- count_df[sub_genes, ]
            spca_obj <- CreateSpatialPCAObject(
                counts = as.matrix(count_sub),
                location = as.matrix(coor_df),
                project = "SpatialPCA",
                gene.type = "custom",
                #sparkversion = "spark",
                #numCores_spark = 5,
                #gene.number = 3000,
                customGenelist = rownames(count_sub),
                min.loctions = 0,
                min.features = 1
            )
            spca_result <- sub_spatialpca(spca_obj, shape, ncs)
            names(spca_result) <- colnames(spca_obj@counts)
            if (!is.null(spca_result)) {
                write.csv(
                    spca_result,
                    paste0(
                        WORKDIR, "results/spatialpca/",
                        sample, "-", method, "-", begin, "_", begin + step,
                        "-", "spatialpca.csv"
                    )
                )
            }

        }
    }

}

read_rank <- function(sample) {
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
        rank_list$spark <- rownames(read_df)
    })

    return(rank_list)
}

# %%
sample <- 151673
shape <- "hexagon"
ncs <- 7
read_files <- read_dlpfc(sample)
sample <- paste0("DLPFC-", sample)
rank_list <- read_rank(sample)
for (step in c(500, 800)) {
    perform_spatialpca(
        read_files$count_df,
        read_files$coor_df,
        rank_list,
        step,
        shape,
        ncs
    )
}

# %%
sample <- "E9.5_E1S1.MOSTA"
shape <- "square"
ncs <- 12
read_files <- read_stereo(sample)
rank_list <- read_rank(sample)
for (step in c(500)) {
    perform_spatialpca(
        read_files$count_df,
        read_files$coor_df,
        rank_list,
        step,
        shape,
        ncs
    )
}

# %%
sample <- "Mouse_brain"
shape <- "square"
ncs <- 20
read_files <- read_stereo(sample)
rank_list <- read_rank(sample)
for (step in c(500)) {
    perform_spatialpca(
        read_files$count_df,
        read_files$coor_df,
        rank_list,
        step,
        shape,
        ncs
    )
}
