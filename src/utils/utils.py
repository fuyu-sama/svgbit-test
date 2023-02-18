from pathlib import Path
import pandas as pd

WORKDIR = Path.joinpath(Path.home(), "workspace/svgbit_local_test")


def replace_label(in_series):
    out_list = []
    names = {}
    j = 0
    for i in in_series.values:
        if i not in names:
            names[i] = j
            j += 1
        out_list.append(names[i])
    return out_list


def read_rank(sample):
    rank_dict = {}

    read_df = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/{sample}/AI.csv"),
        index_col=0,
        header=0,
    ).sort_values(by="AI", ascending=False)
    rank_dict["SVGbit"] = list(read_df.index)

    try:
        read_df = pd.read_csv(
            Path.joinpath(WORKDIR, f"results/somde/{sample}-somde.csv"),
            index_col=0,
            header=0,
        ).sort_values(by="qval")
        rank_dict["SOMDE"] = list(read_df["g"])
    except FileNotFoundError:
        print(f"No SOMDE result for {sample}")
    try:
        read_df = pd.read_csv(
            Path.joinpath(WORKDIR,
                          f"results/spatialde/{sample}-spatialde.csv"),
            index_col=0,
            header=0,
        ).sort_values(by="qval")
        rank_dict["SpatialDE"] = list(read_df["g"])
    except FileNotFoundError:
        print(f"No SpatialDE result for {sample}")
    try:
        read_df = pd.read_csv(
            Path.joinpath(WORKDIR, f"results/spark/{sample}-spark.csv"),
            index_col=0,
            header=0,
        ).sort_values(by="adjusted_pvalue")
        rank_dict["SPARK"] = list(read_df.index)
    except FileNotFoundError:
        print(f"No SPARK result for {sample}")

    return rank_dict
