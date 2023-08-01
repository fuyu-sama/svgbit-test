from pathlib import Path

import pandas as pd
import scanpy as sc

import svgbit as sb

import src.utils.read_label as read_label

WORKDIR = Path.joinpath(Path.home(), "workspace/svgbit_local_test")


def read_mba(gsm, sample):
    # count_df: spots * genes
    count_path = Path.joinpath(
        WORKDIR,
        "data/2020_SciAdv_MouseBrainAtlus/COUNT/",
        f"{gsm}_expr_raw_counts_table_{sample}.tsv.gz",
    )
    count_df = pd.read_csv(
        count_path,
        index_col=0,
        header=0,
        sep="\t",
        compression="gzip",
    ).T
    count_df = count_df[~count_df.index.duplicated(keep="first")].T

    coor_path = Path.joinpath(
        WORKDIR,
        f"data/2020_SciAdv_MouseBrainAtlus/COORDINATES/coordinates_{sample}.tsv/",
    )
    coor_df = pd.read_csv(coor_path, index_col=0, header=0, sep="\t")
    coor_df.columns = ["X", "Y"]

    return (count_df, coor_df)


def read_als(sample):
    # count_df: spots * genes
    count_path = Path.joinpath(
        WORKDIR,
        "data/2019_Science_ALS/COUNT/",
        f"{sample}_stdata_aligned_counts_IDs.txt.gz",
    )
    count_df = pd.read_csv(
        count_path,
        index_col=0,
        header=0,
        sep="\t",
        compression="gzip",
    )
    count_df = count_df[~count_df.index.duplicated(keep="first")].T

    x_coor = []
    y_coor = []
    for i in count_df.index:
        x_coor.append(float(i.split("_")[0]))
        y_coor.append(float(i.split("_")[1]))
    coor_df = pd.DataFrame((x_coor, y_coor)).T
    coor_df.index = count_df.index
    coor_df.columns = [["X", "Y"]]

    return (count_df, coor_df)


def read_stereo(sample):
    # count_df: spots * genes
    adata = sc.read_h5ad(
        f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")
    count_df = pd.DataFrame.sparse.from_spmatrix(
        adata.layers["count"],
        index=adata.obs.index,
        columns=adata.var.index,
    ).sparse.to_dense()

    coor_df = pd.DataFrame(adata.obsm["spatial"])
    coor_df.index = count_df.index
    coor_df.columns = ["X", "Y"]

    return (count_df, coor_df)


def read_10x(sample):
    # count_df: spots * genes
    d = sb.load_10X(f"data/10XGenomics/{sample}/outs")
    count_df = d.count_df.sparse.to_dense()
    coor_df = d.coor_df
    return (count_df, coor_df)


def read_dlpfc(sample):
    # count_df: spots * genes
    count_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-counts.csv",
    )
    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    logcount_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-logcounts.csv",
    )
    logcount_df = pd.read_csv(logcount_path, index_col=0, header=0).T
    coor_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-coor.csv",
    )
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    coor_df.columns = ["X", "Y"]

    label_series = read_label.read_dlpfc_label(sample)
    label_series = label_series.reindex(index=count_df.index)
    na_spots = label_series[label_series.isna()].index
    count_df = count_df.drop(index=na_spots)
    coor_df = coor_df.drop(index=na_spots)
    logcount_df = logcount_df.drop(index=na_spots)

    return (count_df, coor_df, logcount_df)


def read_self(sample):
    path = Path.joinpath(WORKDIR, f"data/mouse-brain-full/{sample}/outs")
    d = sb.load_10X(path, make_sparse=False)
    d = sb.filters.low_variance_filter(d, 0)
    return (d.count_df, d.coordinate_df)
