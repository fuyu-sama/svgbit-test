# svgbit-test

This repository contains test codes used by *The spatiotemporal dynamics of
spatially variable genes in developing mouse brain revealed by a novel
computational scheme*. For detailed information, please refer to our publication.

Paper link:
https://www.nature.com/articles/s41420-023-01569-w

For more information about SVGbit, please visit
https://github.com/CPenglab/svgbit

## File detail

- requirements.txt

    Python packages and version used.

- src/svgbit_pipeline_test

    Codes for run svgbit on several datasets.

- src/utils

    Codes for util functions, including reading functions.

- src/cluster_index.py

    Python code for calculate ARI.

- src/draw_k.py

    Python code for draw correlation between different K values.

- src/run-bayesspace.R

    R code for run BayesSpace clustering with different step length.

- src/run-somde.py
- src/run-spark.R
- src/run-spatialde.py

    Codes for discovering siatially variable genes with SOMDE/SPARK/SpatialDE.

- src/run-spatialpca.R

    R code for run SpatialPCA clustering with different step length. Not used in
    our publication.

- src/test_k.py

    Python code for run SVGbit with different K.
