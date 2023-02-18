from pathlib import Path

import pandas as pd
import scanpy as sc

WORKDIR = Path.joinpath(Path.home(), "workspace/svgbit_local_test")


def strip_zero(str_in):
    str_in = str(str_in)
    while 1:
        if str_in[-1] == "0":
            str_in = str_in[:-1]
        else:
            return str_in


def read_als_label(sample):
    label_path = Path.joinpath(
        WORKDIR,
        "data/2019_Science_ALS/ANNOTATION",
        f"{sample}.tsv.annotations.tsv.gz",
    )
    label_df = pd.read_csv(
        label_path,
        index_col=None,
        header=0,
        sep="\t",
        compression="gzip",
    )
    label_dict = {}
    for index, row in label_df.iterrows():
        spot_name = f"{strip_zero(row['xPos'])}_{strip_zero(row['yPos'])}"
        row = row[4:]
        assert (len(row[row == 1]) == 1)
        label_dict[spot_name] = row[row == 1].index[0]

    return pd.DataFrame.from_dict(label_dict, orient="index")[0]


def read_mba_label(gsm, sample):
    meta_path = Path.joinpath(
        WORKDIR,
        "data/2020_SciAdv_MouseBrainAtlus/META/",
        f"{gsm}_meta_table_{sample}.tsv.gz",
    )
    meta_df = pd.read_csv(
        meta_path,
        index_col=0,
        header=0,
        sep="\t",
        compression="gzip",
    )

    return meta_df["ABA_parent"]


def read_stereo_label(sample):
    adata = sc.read_h5ad(
        f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")

    return adata.obs["annotation"]


def read_dlpfc_label(sample):
    coldata_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-coldata.csv",
    )
    coldata_df = pd.read_csv(coldata_path, index_col=0, header=0)

    return coldata_df["spatialLIBD"]
