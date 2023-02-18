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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import svgbit as sb

import src.utils.read_func as read_func

plt.rcParams.update({"font.size": 18})


def perform_k_coor(count_df, coor_df, sample):
    ai_dict = {}
    ai_df = pd.DataFrame()
    for k in range(4, 13):
        d = sb.STDataset(count_df, coor_df)
        d = sb.filters.low_variance_filter(d)
        d = sb.filters.high_expression_filter(d, 0.99)
        d = sb.normalizers.logcpm_normalizer(d)
        d.acquire_weight(k=k)
        d.acquire_hotspot(cores=5)
        d.acquire_density(cores=5)
        ai_series = d.AI.copy()
        ai_series.name = k
        ai_dict[k] = ai_series
        ai_df = pd.concat([ai_df, ai_series], axis=1)
        del d

    fig, ax = plt.subplots(figsize=(10, 10))
    ai_df.corr(
        method="spearman").to_csv(f"results/k_corr/{sample}-spearman.csv")
    ai_df.corr(method="pearson").to_csv(f"results/k_corr/{sample}-pearson.csv")
    sns.heatmap(
        ai_df.corr(method="spearman"),
        ax=ax,
        annot=True,
        linewidth=.5,
        cmap="Reds",
        vmax=1.1,
        vmin=0.6,
        cbar=False,
    )
    ax.set_xlabel("K")
    ax.set_ylabel("K")
    ax.invert_yaxis()
    fig.savefig(f"results/k_corr/{sample}-spearman.pdf")

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(
        ai_df.corr(method="pearson"),
        ax=ax,
        annot=True,
        linewidth=.5,
        cmap="Reds",
        vmax=1.1,
        vmin=0.6,
        cbar=False,
    )
    ax.set_xlabel("K")
    ax.set_ylabel("K")
    ax.invert_yaxis()
    fig.savefig(f"results/k_corr/{sample}-pearson.pdf")


# %%
sample = 151673
count_df, coor_df, *_ = read_func.read_dlpfc(sample)
perform_k_coor(count_df, coor_df, f"DLPFC-{sample}")
