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
    out_list = pd.Series(out_list, index=in_series.index)
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


def jaccard_index(set_1, set_2):
    return (len(set_1) & len(set_2)) / len(set_1.union(set_2))


def match_labels(labels, results):
    return_series = pd.Series()
    for cluster in set(results.values):
        cluster_spots = set(results[results == cluster].index)
        max_dist, max_dist_region = -1, None
        for region in set(labels.values):
            region_spots = set(labels[labels == region].index)
            dist = jaccard_index(region_spots, cluster_spots)
            if dist > max_dist:
                max_dist, max_dist_region = dist, region
        temp_series = pd.Series(
            [max_dist_region] * len(cluster_spots),
            index=cluster_spots,
        )
        return_series = pd.concat([return_series, temp_series])
    return_series = return_series.reindex(index=results.index)
    return return_series
