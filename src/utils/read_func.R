WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/svgbit_local_test/")
renv::activate(WORKDIR)
reticulate::use_virtualenv(paste0(WORKDIR, "venv"))

read_als <- function(sample) {
    # count_df: genes * spots
    count_path <- paste0(
        WORKDIR,
        "data/2019_Science_ALS/COUNT/", sample, "_stdata_aligned_counts_IDs.txt.gz"
    )
    count_df <- read.csv(count_path, sep = "\t", row.names=NULL)
    gene_names <- count_df$row.names
    count_df <- count_df[, -1]

    x_coor <- c()
    y_coor <- c()
    for (spot_name in colnames(count_df)) {
        spot_name <- strsplit(spot_name, "_")[[1]]
        x_coor <- c(x_coor, as.numeric(strsplit(spot_name[1], "X")[[1]][2]))
        y_coor <- c(y_coor, as.numeric(spot_name[2]))
    }
    coor_df <- data.frame(X = x_coor, Y = y_coor)
    rownames(coor_df) <- colnames(count_df)

    return(list(count_df = count_df, coor_df = coor_df, gene_names = gene_names))
}

read_mba <- function(gsa, sample) {
    # count_df: genes * spots
    count_path <- paste0(
        WORKDIR,
        "data/2020_SciAdv_MouseBrainAtlus/COUNT/",
        gsm, "_expr_raw_counts_table_", sample, ".tsv.gz"
    )
    count_df <- read.csv(count_path, sep = "\t", row.names = 1)
    count_df <- t(count_df)

    x_coor <- c()
    y_coor <- c()
    for (spot_name in colnames(count_df)) {
        spot_name <- strsplit(spot_name, "x")[[1]]
        x_coor <- c(x_coor, as.numeric(strsplit(spot_name[1], "_")[[1]][2]))
        y_coor <- c(y_coor, as.numeric(spot_name[2]))
    }
    coor_df <- data.frame(X = x_coor, Y = y_coor)
    rownames(coor_df) <- colnames(count_df)

    return(list(count_df = count_df, coor_df = coor_df))
}

read_stereo <- function(sample) {
    # count_df: genes * spots
    read_path <- paste0(
        WORKDIR,
        "data/2022_Cell_Stereo-seq/STDS0000058/stomics/", sample, ".h5ad"
    )
    adata <- anndata::read_h5ad(read_path)
    count_df <- t(as.data.frame(as.matrix(adata$layers["count"])))
    logcounts_df <- t(as.data.frame(as.matrix(adata$X)))
    coor_df <- as.data.frame(adata$obsm$spatial)
    colnames(coor_df) <- c("X", "Y")
    rownames(coor_df) <- colnames(count_df)

    return(list(
            count_df = count_df, 
            logcounts_df = logcounts_df, 
            coor_df = coor_df))
}

read_dlpfc <- function(sample) {
    # count_df: genes * spots
    count_path <- paste0(
        WORKDIR, "data/spatialLIBD/", sample, "-counts.csv"
    )
    logcounts_path <- paste0(
        WORKDIR, "data/spatialLIBD/", sample, "-logcounts.csv"
    )
    coor_path <- paste0(
        WORKDIR, "data/spatialLIBD/", sample, "-coor.csv"
    )
    coldata_path <- paste0(
        WORKDIR, "data/spatialLIBD/", sample, "-coldata.csv"
    )

    count_df <- read.csv(count_path, row.names = 1, check.names = FALSE)
    logcounts_df <- read.csv(logcounts_path, row.names = 1, check.names = FALSE)
    coor_df <- read.csv(coor_path, row.names = 1)
    colnames(coor_df) <- c("X", "Y")

    coldata <- read.csv(coldata_path, row.names = 1)
    coldata <- coldata[complete.cases(coldata[, "spatialLIBD"]),]
    count_df <- count_df[, rownames(coldata)]
    logcounts_df <- logcounts_df[, rownames(coldata)]
    coor_df <- coor_df[rownames(coldata), ]
    array_df <- coldata[, c("array_row", "array_col")]
    colnames(array_df) <- c("row", "col")

    return(list(
            count_df = count_df,
            logcounts_df = logcounts_df,
            coor_df = coor_df,
            array_df = array_df))
}
