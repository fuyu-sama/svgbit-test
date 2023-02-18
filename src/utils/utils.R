split_coor <- function(spot_names) {
    coor_df <- lapply(spot_names, function(spot_name) {
        return(as.numeric(strsplit(spot_name, "_")[[1]]))
    })
    coor_df <- t(as.data.frame(coor_df))
    colnames(coor_df) <- c("X", "Y")
    rownames(coor_df) <- spot_names
    return(coor_df)
}
