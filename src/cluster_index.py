#! venv/bin/python
# -*- encoding: utf-8 -*-

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#  Codes are far away from bugs with Buddha's bless
#

# %%
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

import src.utils.read_func as read_func
import src.utils.read_label as read_label
import src.utils.read_image as read_image
from src.utils.utils import WORKDIR, replace_label, read_rank

plt.rcParams.update({"font.size": 18})

score_funcs = [
    ("V-measure", metrics.v_measure_score),
    ("Rand_index", metrics.rand_score),
    ("ARI", metrics.adjusted_rand_score),
    ("MI", metrics.mutual_info_score),
    ("NMI", metrics.normalized_mutual_info_score),
    ("AMI", metrics.adjusted_mutual_info_score),
    ("FMI", metrics.fowlkes_mallows_score),
]

seed = 42


def perform_bayesspace(sample, label_series, step):
    bayesspace_result = {}

    for method in ["SVGbit", "SOMDE", "SpatialDE", "spark", "SPARK"]:
        matrix_dict = {"ncs": [], "begin": [], "end": []}
        for i in score_funcs:
            matrix_dict[i[0]] = []
        cluster_results = {}
        flag = 0
        for begin in range(0, step * 5, step):
            read_path = Path.joinpath(
                WORKDIR,
                "results/bayesspace/",
                f"seed-{seed}/",
                f"{sample}-{method}-{begin}_{begin + step}-bayesspace.csv",
            )
            try:
                read_df = pd.read_csv(read_path, index_col=0, header=0)
                cluster_result = read_df["spatial.cluster"]
                flag = 1
            except FileNotFoundError:
                print(
                    "No BayesSpace result for",
                    f"{method}-{begin}_{begin + step}",
                )
                continue
            cluster_results[f"{begin} - {begin + step}"] = cluster_result
            matrix_dict["ncs"].append(len(set(cluster_result)))
            matrix_dict["begin"].append(begin)
            matrix_dict["end"].append(begin + step)
            label_series = label_series.reindex(index=cluster_result.index)
            for i in score_funcs:
                matrix_dict[i[0]].append(i[1](label_series, cluster_result))
        if flag:
            if method == "spark":
                method = "SPARK"
            bayesspace_result[method] = {
                "matrix_df": pd.DataFrame.from_dict(matrix_dict),
                "result": cluster_results,
            }

    return bayesspace_result


def perform_spatialpca(sample, label_series, step):
    spca_result = {}

    for method in ["SVGbit", "SOMDE", "SpatialDE", "SPARK"]:
        matrix_dict = {"ncs": [], "begin": [], "end": []}
        for i in score_funcs:
            matrix_dict[i[0]] = []
        cluster_results = {}
        flag = 0
        for begin in range(0, step * 5, step):
            read_path = Path.joinpath(
                WORKDIR,
                "results/spatialpca/",
                f"{sample}-{method}-{begin}_{begin + step}-spatialpca.csv",
            )
            try:
                read_df = pd.read_csv(read_path, index_col=0, header=0)
                cluster_result = read_df["x"]
                flag = 1
            except FileNotFoundError:
                print(
                    "No SpatialPCA result for",
                    f"{method}-{begin}_{begin + step}",
                )
                continue
            cluster_results[f"{begin} - {begin + step}"] = cluster_result
            matrix_dict["ncs"].append(len(set(cluster_result)))
            matrix_dict["begin"].append(begin)
            matrix_dict["end"].append(begin + step)
            label_series = label_series.reindex(index=cluster_result.index)
            for i in score_funcs:
                matrix_dict[i[0]].append(i[1](label_series, cluster_result))
        if flag:
            spca_result[method] = {
                "matrix_df": pd.DataFrame.from_dict(matrix_dict),
                "result": cluster_results,
            }

    return spca_result


def _perform_kmeans(rank_list, count_df, label_series, step):
    result_dict = {"ncs": [], "begin": [], "end": []}
    kmeans_results = {}
    for i in score_funcs:
        result_dict[i[0]] = []
    ncs = len(set(label_series))
    for begin in range(0, 2000 + 1, step):
        count_sub = count_df.reindex(columns=rank_list[begin:begin + step])
        cluster_result = KMeans(n_clusters=ncs).fit_predict(
            PCA(n_components=40).fit_transform(count_sub))
        cluster_result = pd.Series(cluster_result, index=count_df.index)
        kmeans_results[f"{begin} - {begin + step}"] = cluster_result
        result_dict["ncs"].append(ncs)
        result_dict["begin"].append(begin)
        result_dict["end"].append(begin + step)
        for i in score_funcs:
            result_dict[i[0]].append(i[1](label_series, cluster_result))

    return (pd.DataFrame.from_dict(result_dict), kmeans_results)


def perform_kmeans(rank_dict, *args, **kwargs):
    kmeans_result_dict = {}
    for method in rank_dict:
        matrix_df, kmeans_result = _perform_kmeans(
            rank_dict[method],
            *args,
            **kwargs,
        )
        kmeans_result_dict[method] = {
            "matrix_df": matrix_df,
            "result": kmeans_result
        }
    return kmeans_result_dict


def draw_matrix(result_dict, sample, step, title, matrix="ARI"):
    draw_df = pd.DataFrame()
    for method in result_dict:
        draw_sub = pd.Series()
        result_df = result_dict[method]["matrix_df"]
        for begin in range(0, step * 5, step):
            line = result_df[result_df["begin"] == begin]
            line = line[line["end"] == begin + step].iloc[0:1, ][matrix]
            line.index = [f"{begin + 1} - {begin + step}"]
            draw_sub = pd.concat([draw_sub, line])
            draw_sub.name = method
        draw_df = pd.concat([draw_df, draw_sub], axis=1)

    fig, ax = plt.subplots(figsize=(10, 10))
    plot = ax.plot(draw_df, marker="o")
    plt.setp(ax.get_xticklabels(), rotation=45)
    ax.legend(plot, draw_df.columns)
    ax.set_ylim([0, 0.6])
    fig.savefig(
        f"results/cluster_index/{sample}-step{step}-{title}-{matrix}.svg",
        bbox_inches="tight",
    )
    plt.close(fig)


def draw_cluster(result_dict, label_series, s, title, he_image=None):
    for method in result_dict:
        ncol = len(result_dict[method]["result"]) + 1
        fig, axes = plt.subplots(1, ncol, figsize=(ncol * 10, 10))
        [ax.axis("off") for ax in axes]
        if he_image is not None:
            [ax.imshow(he_image) for ax in axes]
        axes[0].scatter(
            coor_df["X"],
            coor_df["Y"],
            s=s,
            cmap="tab20",
            c=replace_label(label_series),
        )
        axes[0].set_title("Annotation region")
        for ax, rank in zip(axes[1:], result_dict[method]["result"]):
            ax.scatter(
                coor_df["X"],
                coor_df["Y"],
                s=s,
                cmap="tab20",
                c=result_dict[method]["result"][rank],
            )
            ax.set_title(rank)
        fig.savefig(
            f"results/cluster_index/{sample}-step{step}-{method}-{title}.svg",
            bbox_inches="tight",
        )
        plt.close(fig)


def draw_cluster_separate(
    result_dict,
    label_series,
    s,
    title,
    he_image=None,
    im_format="jpg",
):
    for method in result_dict:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.axis("off")
        ax.imshow(he_image) if he_image is not None else None
        ax.scatter(
            coor_df["X"],
            coor_df["Y"],
            s=s,
            cmap="tab20",
            c=replace_label(label_series),
        )
        ax.set_title("Annotation region")
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/cluster_index/",
                f"{sample}-step{step}-{method}-{title}-all.{im_format}",
            ),
            bbox_inches="tight",
        )
        plt.close(fig)
        for rank in result_dict[method]["result"]:
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(he_image) if he_image is not None else None
            ax.axis("off")
            ax.scatter(
                coor_df["X"],
                coor_df["Y"],
                s=s,
                cmap="tab20",
                c=result_dict[method]["result"][rank],
            )
            ax.set_title(rank)
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster_index/",
                    f"{sample}-step{step}-{method}-{title}-{rank}.{im_format}",
                ),
                bbox_inches="tight",
            )
            plt.close(fig)


def draw_tsne(rank_dict, count_df, label_series, step):
    label_list = replace_label(label_series)
    full_tsne_result = TSNE().fit_transform(
        PCA(n_components=40).fit_transform(count_df))
    for method, rank_list in rank_dict.items():
        ranges = range(0, step * 5, step)
        ncol = len(ranges)
        fig, axes = plt.subplots(1, ncol, figsize=(ncol * 10, 10))
        axes[0].scatter(
            full_tsne_result[:, 0],
            full_tsne_result[:, 1],
            c=label_list,
            cmap="tab20",
        )
        axes[0].set_title("All genes")
        for begin, ax in zip(ranges, axes[1:]):
            count_sub = count_df.reindex(columns=rank_list[begin:begin + step])
            tsne_result = TSNE().fit_transform(count_sub)
            ax.set_title(f"{begin + 1} - {begin + step}")
            ax.scatter(
                tsne_result[:, 0],
                tsne_result[:, 1],
                c=label_list,
                cmap="tab20",
            )
        fig.savefig(
            f"results/cluster_index/{sample}-step{step}-{method}-pcatsne.pdf",
            bbox_inches="tight",
        )
        plt.close(fig)


def draw_umap(rank_dict, count_df, label_series, step):
    label_list = replace_label(label_series)
    full_umap_result = UMAP().fit_transform(
        PCA(n_components=40).fit_transform(count_df))
    for method, rank_list in rank_dict.items():
        ranges = range(0, step * 5, step)
        ncol = len(ranges)
        fig, axes = plt.subplots(1, ncol, figsize=(ncol * 10, 10))
        axes[0].scatter(
            full_umap_result[:, 0],
            full_umap_result[:, 1],
            c=label_list,
            cmap="tab20",
        )
        axes[0].set_title("All genes")
        for begin, ax in zip(ranges, axes[1:]):
            count_sub = count_df.reindex(columns=rank_list[begin:begin + step])
            tsne_result = UMAP().fit_transform(count_sub)
            ax.set_title(f"{begin + 1} - {begin + step}")
            ax.scatter(
                tsne_result[:, 0],
                tsne_result[:, 1],
                c=label_list,
                cmap="tab20",
            )
        fig.savefig(
            f"results/cluster_index/{sample}-step{step}-{method}-pcaumap.pdf",
            bbox_inches="tight",
        )
        plt.close(fig)


# %%
sample = "E9.5_E1S1.MOSTA"
count_df, coor_df = read_func.read_stereo(sample)
label_series = read_label.read_stereo_label(sample)
coor_df = coor_df.reindex(index=count_df.index)
label_series = label_series.reindex(index=coor_df.index)
rank_dict = read_rank(sample)

for step in (500, ):
    bayesspace_result = perform_bayesspace(sample, label_series, step)
    draw_matrix(bayesspace_result, sample, step, "BayesSpace", "ARI")
    draw_cluster(bayesspace_result, label_series, 64, "BayesSpace")
    draw_tsne(rank_dict, count_df, label_series, step)

# %%
sample = 151673
count_df, coor_df, *_ = read_func.read_dlpfc(sample)
label_series = read_label.read_dlpfc_label(sample).reindex(
    index=count_df.index)
coor_df = coor_df.reindex(index=count_df.index)
he_image = read_image.read_dlpfc_image(sample)
sample = f"DLPFC-{sample}"
rank_dict = read_rank(sample)

for step in (500, ):
    bayesspace_result = perform_bayesspace(sample, label_series, step)
    draw_matrix(bayesspace_result, sample, step, "BayesSpace", "ARI")
    draw_cluster(bayesspace_result, label_series, 16, "BayesSpace", he_image)
    kmeans_result = perform_kmeans(rank_dict, count_df, label_series, step)
    draw_cluster(kmeans_result, label_series, 16, "kmeans", he_image)

# %%
sample = "Mouse_brain"
count_df, coor_df = read_func.read_stereo(sample)
label_series = read_label.read_stereo_label(sample)
coor_df = coor_df.reindex(index=count_df.index)
label_series = label_series.reindex(index=coor_df.index)
rank_dict = read_rank(sample)

for step in (500, ):
    bayesspace_result = perform_bayesspace(sample, label_series, step)
    draw_matrix(bayesspace_result, sample, step, "BayesSpace", "ARI")
    draw_cluster_separate(
        {"SVGbit": bayesspace_result["SVGbit"]},
        label_series,
        16,
        "BayesSpace",
        im_format="png",
    )
