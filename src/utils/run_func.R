run_spark <- function(count_df, coor_df) {
    # count_df: genes * spots
    library(SPARK)
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
