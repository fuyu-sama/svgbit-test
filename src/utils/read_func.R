WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "svgbit-test")
reticulate::use_virtualenv(fs::path(WORKDIR, "venv"))

read_mba <- function(sample) {
    # count_df: genes * spots
    gsm <- mba_files[[sample]]
    count_path <- fs::path(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "COUNT",
        paste0(gsm, "_expr_raw_counts_table_", sample, ".tsv.gz")
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
    coor_df <- data.frame("X" = x_coor, "Y" = y_coor)
    array_df <- data.frame("array_col" = x_coor, "array_row" = y_coor)
    rownames(coor_df) <- colnames(count_df)

    returns <- list(
        "sample" = sample,
        "source" = "mba",
        "count_df" = count_df,
        "coor_df" = coor_df,
        "array_df" = array_df,
        "label" = read_mba_label(sample)
    )
    return(returns)
}

read_mba_label <- function(sample) {
    meta_path <- fs::path(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "META",
        paste0(mba_files[[sample]], "_meta_table_", sample, ".tsv.gz")
    )

    meta_df <- read.table(
        meta_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
    )

    returns <- meta_df$ABA_parent
    names(returns) <- rownames(returns)
    return(returns)
}

list_mba_files <- function() {
    files <- fs::dir_ls(
        fs::path(WORKDIR, "data", "2020_SciAdv_MouseBrainAtlus", "COUNT")
    )
    result_dict <- list()

    for (file in files) {
        file_name <- fs::path_file(file)

        if (startsWith(file_name, "GSM")) {
            match <- regexec("^(GSM\\d+)_.*_(\\d{2}A)\\.tsv\\.gz$", file_name)
            match_result <- regmatches(file_name, match)

            if (length(match_result[[1]]) > 0) {
                gsm <- match_result[[1]][2]
                suffix <- match_result[[1]][3]
                result_dict[[suffix]] <- gsm
            }
        }
    }
    return(result_dict)
}
mba_files <- list_mba_files()

read_stereo <- function(sample) {
    # count_df: genes * spots
    read_path <- fs::path(
        WORKDIR,
        "data",
        "2022_Cell_Stereo-seq",
        "STDS0000058",
        "stomics",
        paste0(sample, ".h5ad")
    )
    adata <- anndata::read_h5ad(read_path)
    count_df <- t(as.data.frame(as.matrix(adata$layers["count"])))
    logcount_df <- t(as.data.frame(as.matrix(adata$X)))
    coor_df <- as.data.frame(adata$obsm$spatial)
    colnames(coor_df) <- c("X", "Y")
    rownames(coor_df) <- colnames(count_df)
    array_df <- coor_df
    colnames(array_df) <- c("array_col", "array_row")

    returns <- list(
        "sample" = sample,
        "source" = "stereo",
        "count_df" = count_df,
        "logcount_df" = logcount_df,
        "coor_df" = coor_df,
        "array_df" = array_df
    )
    return(returns)
}

read_dlpfc <- function(sample) {
    # count_df: genes * spots
    count_path <- fs::path(
        WORKDIR, "data", "spatialLIBD", paste0(sample, "-counts.csv")
    )
    logcount_path <- fs::path(
        WORKDIR, "data", "spatialLIBD", paste0(sample, "-logcounts.csv")
    )
    coor_path <- fs::path(
        WORKDIR, "data", "spatialLIBD", paste0(sample, "-coor.csv")
    )
    coldata_path <- fs::path(
        WORKDIR, "data", "spatialLIBD", paste0(sample, "-coldata.csv")
    )

    count_df <- read.csv(count_path, row.names = 1, check.names = FALSE)
    logcount_df <- read.csv(logcount_path, row.names = 1, check.names = FALSE)
    coor_df <- read.csv(coor_path, row.names = 1)
    colnames(coor_df) <- c("X", "Y")

    coldata <- read.csv(coldata_path, row.names = 1)
    coldata <- coldata[complete.cases(coldata[, "spatialLIBD"]),]
    count_df <- count_df[, rownames(coldata)]
    logcount_df <- logcount_df[, rownames(coldata)]
    coor_df <- coor_df[rownames(coldata), ]
    array_df <- coldata[, c("array_row", "array_col")]

    returns <- list(
        "sample" = sample,
        "source" = "dlpfc",
        "count_df" = count_df,
        "logcount_df" = logcount_df,
        "coor_df" = coor_df,
        "array_df" = array_df
    )
    return(returns)
}

read_iimpact <- function(sample) {
    # count_df: genes * spots
    sample_dir <- paste0(gsub("_", " ", sample), " FFPE data")
    count_path <- fs::path(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        paste0("10x_", sample, "_ffpe_count.csv")
    )
    count_df <- read.csv(count_path, row.names = 1, check.names = FALSE)
    count_df <- count_df[, colSums(count_df) > 0]
    coor_path <- fs::path(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        paste0("10x_", sample, "_ffpe_loc.csv")
    )
    coor_df <- read.csv(coor_path, row.names = 1, check.names = FALSE)
    coor_df$index <- paste(coor_df$x, coor_df$y, sep = "_")
    rownames(coor_df) <- coor_df$index
    rownames(count_df) <- coor_df$index
    coor_df <- coor_df[, c("x", "y")]
    colnames(coor_df) <- c("X", "Y")

    label <- read_iimpact_label(sample)
    count_df <- count_df[names(label), ]
    coor_df <- coor_df[names(label), ]
    count_df <- t(count_df)

    dist_matrix <- as.matrix(dist(coor_df))
    min_distance <- min(dist_matrix[dist_matrix > 0])
    scaling_factor <- 1 / min_distance
    array_df <- ceiling(coor_df * scaling_factor)
    min_coords <- apply(array_df, 2, min)
    array_df <- as.data.frame(array_df - min_coords + c(1, 1))
    colnames(array_df) <- c("array_col", "array_row")
    rownames(array_df) <- rownames(coor_df)

    returns <- list(
        "sample" = sample,
        "source" = "iimpact",
        "count_df" = count_df,
        "coor_df" = coor_df,
        "array_df" = array_df,
        "label" = label
    )
    return(returns)
}

read_iimpact_label <- function(sample) {
    sample_dir <- paste0(gsub("_", " ", sample), " FFPE data")
    coor_path <- fs::path(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        paste0("10x_", sample, "_ffpe_manual_annotation.csv")
    )
    coor_df <- read.csv(coor_path, row.names = 1, check.names = FALSE)
    coor_df$index <- paste(coor_df$x, coor_df$y, sep = "_")
    rownames(coor_df) <- coor_df$index
    returns <- coor_df$annotation
    names(returns) <- coor_df$index
    return(returns)
}

read_self <- function(sample) {
    # count_df: genes * spots
    file_name <- "filtered_feature_bc_matrix"
    count_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", sample, "outs", file_name, "matrix.mtx.gz"
    )
    features_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", sample, "outs", file_name, "features.tsv.gz"
    )
    barcodes_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", sample, "outs", file_name, "barcodes.tsv.gz"
    )
    coor_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", sample, "outs",
        "spatial", "tissue_positions_list.csv"
    )

    count_df <- as.matrix(Matrix::readMM(gzfile(count_path)))
    features <- read.table(gzfile(features_path), header = FALSE, sep = "\t")$V2
    barcodes <- read.table(gzfile(barcodes_path), header = FALSE, sep = "\t")$V1
    rownames(count_df) <- make.unique(features)
    colnames(count_df) <- barcodes

    coor_df <- read.csv(coor_path, row.names = 1, header = FALSE)
    coor_df <- coor_df[colnames(count_df), ]
    array_df <- coor_df[, c("V3", "V2")]
    coor_df <- coor_df[, c("V5", "V4")]
    colnames(coor_df) <- c("X", "Y")
    colnames(array_df) <- c("array_col", "array_row")

    returns <- list(
        "sample" = sample,
        "source" = "self",
        "count_df" = count_df,
        "coor_df" = coor_df,
        "array_df" = array_df
    )
    return(returns)
}

read_self.2 <- function(sample) {
    # count_df: genes * spots
    count_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", "Data", "scale_df", "raw",
        paste0(sample, "-raw.csv")
    )
    coor_path <- fs::path(
        WORKDIR, "data", "mouse-brain-full", "Data", "coor_df",
        fs::path(sample, "-coor.csv")
    )

    count_df <- read.csv(count_path, row.names = 1, check.names = FALSE)
    coor_df <- read.csv(coor_path, row.names = 1)
    colnames(coor_df) <- c("X", "Y")
    coor_df <- coor_df[colnames(count_df), ]

    returns <- list(
        "sample" = sample,
        "source" = "self",
        "count_df" = count_df,
        "coor_df" = coor_df
    )
    return(returns)
}

read_gbm <- function(sample) {
    # count_df: genes * spots
    file_name <- "filtered_feature_bc_matrix"
    if (sample == "21B-603-5") file_name <- "raw_feature_bc_matrix"
    count_path <- fs::path(
        WORKDIR, "data", "gbm", sample, "outs", file_name, "matrix.mtx.gz"
    )
    features_path <- fs::path(
        WORKDIR, "data", "gbm", sample, "outs", file_name, "features.tsv.gz"
    )
    barcodes_path <- fs::path(
        WORKDIR, "data", "gbm", sample, "outs", file_name, "barcodes.tsv.gz"
    )
    coor_path <- fs::path(
        WORKDIR, "data", "gbm", sample, "outs", "spatial", "tissue_positions.csv"
    )

    count_df <- as.matrix(Matrix::readMM(gzfile(count_path)))
    features <- read.table(gzfile(features_path), header = FALSE, sep = "\t")$V2
    barcodes <- read.table(gzfile(barcodes_path), header = FALSE, sep = "\t")$V1
    rownames(count_df) <- make.unique(features)
    colnames(count_df) <- barcodes

    coor_df <- read.csv(coor_path, row.names = 1, header = TRUE)
    coor_df <- coor_df[colnames(count_df), ]
    array_df <- coor_df[, c("array_col", "array_row")]
    coor_df <- coor_df[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
    colnames(coor_df) <- c("X", "Y")

    returns <- list(
        "sample" = sample,
        "source" = "gbm",
        "count_df" = count_df,
        "coor_df" = coor_df,
        "array_df" = array_df
    )
    return(returns)
}
