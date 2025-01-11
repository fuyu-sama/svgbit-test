run_spark <- function(count_df, coor_df) {
    # count_df: genes * spots
    library('SPARK')
    spark_obj <- CreateSPARKObject(
        counts = count_df,
        location = coor_df,
    )
    spark_obj@lib_size <- apply(spark_obj@counts, 2, sum)
    spark_obj <- spark.vc(
        spark_obj,
        covariates = NULL,
        lib_size = spark_obj@lib_size,
        num_core = 5,
        verbose = FALSE
    )
    spark_obj <- spark.test(
        spark_obj,
        check_positive = TRUE,
        verbose = FALSE
    )
    spark_result <- subset(
        spark_obj@res_mtest,
        select = c("combined_pvalue", "adjusted_pvalue")
    )
    spark_result <- spark_result[order(spark_result, decreasing = FALSE), ]
    spark_result <- na.omit(spark_result)

    return(spark_result)
}

source(fs::path(WORKDIR, "src", "utils", "heartsvg_fast.R"))
run_heartsvg <- function(count_df, coor_df) {
    # count_df: genes * spots
    library('HEARTSVG')
    heartsvg.input <- cbind(coor_df[c(2, 1)], t(count_df))
    colnames(heartsvg.input)[1:2] <- c("row", "col")
    result <- heartsvg_fast(heartsvg.input)
    return(result)
}

run_despace <- function(count_df, coor_df, platform = "Visium", ncs = 10) {
    library('DESpace')
    sce <- SingleCellExperiment(
        assays = list(counts = as(as.matrix(count_df), "dgCMatrix")),
        colData = coor_df
    )
    sce <- spatialPreprocess(
        sce,
        platform = platform,
        n.PCs = 15,
        n.HVGs = 2000,
        skip.PCA = FALSE,
        log.normalize = TRUE,
        assay.type = "logcounts"
    )
    sce <- spatialCluster(
        sce,
        q = ncs,
        platform = platform,
        init.method = "mclust",
        model = "t"
    )
    results <- DESpace_test(sce, spatial_cluster = "spatial.cluster", verbose = FALSE)
    return(results$gene_results)
}

run_meringue <- function(count_df, coor_df) {
    # count_df: genes * spots
    counts <- cleanCounts(
        counts = as.matrix(count_df),
        min.reads = 100,
        min.lib.size = 100,
        plot = FALSE,
        verbose = FALSE
    )
    pos <- coor_df[colnames(counts), ]
    mat <- normalizeCounts(counts = counts, log = FALSE, verbose = FALSE)
    knn_result <- FNN::knn.dist(pos, k = 8)
    knn_thresholds <- mean(knn_result[, 1])
    w <- getSpatialNeighbors(pos, filterDist = knn_thresholds)
    I <- getSpatialPatterns(mat, w, verbose = FALSE)
    return(I)
}
