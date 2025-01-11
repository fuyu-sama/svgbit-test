import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional
from PIL import Image

import pandas as pd
import scanpy as sc

import svgbit as sb

Image.MAX_IMAGE_PIXELS = None
WORKDIR = Path.joinpath(Path.home(), "workspace", "svgbit-test")


@dataclass
class ReadData:
    sample: str
    source: str
    count_df: pd.DataFrame
    coor_df: pd.DataFrame
    array_df: pd.DataFrame = None
    logcount_df: Optional[pd.DataFrame] = None
    label: Optional[pd.Series] = None
    image: Optional[Callable[[str], Image.Image]] = None


def read_mba(sample):
    # count_df: spots * genes
    count_path = Path.joinpath(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "COUNT",
        f"{mba_files[sample]}_expr_raw_counts_table_{sample}.tsv.gz",
    )
    count_df = pd.read_csv(
        count_path,
        index_col=0,
        header=0,
        sep="\t",
        compression="gzip",
    ).T
    count_df = count_df[~count_df.index.duplicated(keep="first")].T
    count_df = count_df.loc[:, count_df.sum() > 0]

    coor_path = Path.joinpath(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "COORDINATES",
        f"coordinates_{sample}.tsv",
    )
    coor_df = pd.read_csv(coor_path, index_col=0, header=0, sep="\t")
    coor_df.columns = ["X", "Y"]

    returns = ReadData(
        sample=sample,
        source="mba",
        count_df=count_df,
        coor_df=coor_df,
        label=read_mba_label(sample),
        image=read_mba_image,
    )

    return returns


def read_mba_label(sample):
    meta_path = Path.joinpath(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "META",
        f"{mba_files[sample]}_meta_table_{sample}.tsv.gz",
    )
    meta_df = pd.read_csv(
        meta_path,
        index_col=0,
        header=0,
        sep="\t",
        compression="gzip",
    )
    returns = meta_df["ABA_parent"]
    returns.name = "label"

    return returns


def read_mba_image(gsm, sample):
    image_path = Path.joinpath(
        WORKDIR, f"data/2020_SciAdv_MouseBrainAtlus/HE/{gsm}_HE_{sample}.jpg")
    return Image.open(image_path)


def list_mba_files():
    folder_path = Path.joinpath(
        WORKDIR,
        "data",
        "2020_SciAdv_MouseBrainAtlus",
        "COUNT",
    )
    result_dict = {}
    for file in folder_path.iterdir():
        if file.is_file() and file.name.startswith("GSM"):
            match = re.match(r"(GSM\d+)_.*_(\d{2}A)\.tsv\.gz", file.name)
            if match:
                gsm = match.group(1)
                suffix = match.group(2)
                result_dict[suffix] = gsm
    return result_dict


mba_files = list_mba_files()


def read_stereo(sample):
    # count_df: spots * genes
    adata = sc.read_h5ad(
        f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")
    count_df = pd.DataFrame.sparse.from_spmatrix(
        adata.layers["count"],
        index=adata.obs.index,
        columns=adata.var.index,
    ).sparse.to_dense()
    count_df = count_df.loc[:, count_df.sum() > 0]

    coor_df = pd.DataFrame(adata.obsm["spatial"])
    coor_df.index = count_df.index
    coor_df.columns = ["X", "Y"]

    returns = ReadData(
        sample=sample,
        source="stereo",
        count_df=count_df,
        coor_df=coor_df,
        label=read_stereo_label(sample),
    )

    return returns


def read_stereo_label(sample):
    adata = sc.read_h5ad(
        f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")

    returns = adata.obs["annotation"]
    returns.name = "label"

    return returns


def read_dlpfc(sample):
    # count_df: spots * genes
    count_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-counts.csv",
    )
    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    count_df = count_df.loc[:, count_df.sum() > 0]
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

    label_series = read_dlpfc_label(sample)
    label_series = label_series.reindex(index=count_df.index)
    na_spots = label_series[label_series.isna()].index
    count_df = count_df.drop(index=na_spots)
    coor_df = coor_df.drop(index=na_spots)
    logcount_df = logcount_df.drop(index=na_spots)

    returns = ReadData(
        sample=sample,
        source="dlpfc",
        count_df=count_df,
        coor_df=coor_df,
        logcount_df=logcount_df,
        label=read_dlpfc_label(sample),
        image=read_dlpfc_image,
    )

    return returns


def read_dlpfc_label(sample):
    coldata_path = Path.joinpath(
        WORKDIR,
        f"data/spatialLIBD/{sample}-coldata.csv",
    )
    coldata_df = pd.read_csv(coldata_path, index_col=0, header=0)

    returns = coldata_df["spatialLIBD"]
    returns.name = "labe"

    return returns


def read_dlpfc_image(sample):
    image_path = Path.joinpath(WORKDIR,
                               f"data/spatialLIBD/{sample}_full_image.tif")
    return Image.open(image_path)


def read_iimpact(sample):
    # count_df: spots * genes
    sample_dir = f"{sample.replace('_', ' ')} FFPE data"
    count_path = Path.joinpath(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        f"10x_{sample}_ffpe_count.csv",
    )
    count_df = pd.read_csv(count_path, index_col=0, header=0)
    count_df = count_df.loc[:, count_df.sum() > 0]
    coor_path = Path.joinpath(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        f"10x_{sample}_ffpe_loc.csv",
    )
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    coor_df.index = [f"{data.x}_{data.y}" for _, data in coor_df.iterrows()]
    coor_df.columns = ["X", "Y"]
    count_df.index = coor_df.index

    label = read_iimpact_label(sample)
    count_df = count_df.reindex(index=label.index)
    coor_df = coor_df.reindex(index=label.index)

    returns = ReadData(
        sample=sample,
        source="iimpact",
        count_df=count_df,
        coor_df=coor_df,
        label=label,
    )

    return returns


def read_iimpact_label(sample):
    sample_dir = f"{sample.replace('_', ' ')} FFPE data"
    coor_path = Path.joinpath(
        WORKDIR,
        "data",
        "2024_GB_iIMPACT",
        sample_dir,
        f"10x_{sample}_ffpe_manual_annotation.csv",
    )
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    coor_df.index = [f"{data.x}_{data.y}" for _, data in coor_df.iterrows()]

    returns = coor_df["annotation"]
    returns.name = "label"

    return returns


def read_10x(sample):
    # count_df: spots * genes
    d = sb.load_10X(f"data/10XGenomics/{sample}/outs")
    count_df = d.count_df.sparse.to_dense()
    count_df = count_df.loc[:, count_df.sum() > 0]
    coor_df = d.coor_df
    returns = ReadData(
        sample=sample,
        source="10X",
        count_df=count_df,
        coor_df=coor_df,
    )
    return returns


def read_self(sample):
    path = Path.joinpath(WORKDIR, f"data/mouse-brain-full/{sample}/outs")
    d = sb.load_10X(path, make_sparse=False)
    d = sb.filters.low_variance_filter(d, 0)
    returns = ReadData(
        sample=sample,
        source="self",
        count_df=d.count_df,
        coor_df=d.coordinate_df,
    )
    return returns


def read_gbm(sample):
    path = Path.joinpath(WORKDIR, "data", "gbm", sample, "outs")
    if sample == "21B-603-5":
        file_name = "raw_feature_bc_matrix"
    else:
        file_name = "filtered_feature_bc_matrix"
    d = sb.load_10X(path, file_name=file_name, make_sparse=False)
    d = sb.filters.low_variance_filter(d, 0)
    returns = ReadData(
        sample=sample,
        source="self",
        count_df=d.count_df,
        coor_df=d.coordinate_df,
    )
    return returns


def strip_zero(str_in):
    str_in = str(str_in)
    while 1:
        if str_in[-1] == "0":
            str_in = str_in[:-1]
        else:
            return str_in
