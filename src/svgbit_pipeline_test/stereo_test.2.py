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
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

import svgbit as sb

plt.rcParams.update({"font.size": 18})

# %%
d = sb.load_anndata_h5(
    "data/2022_Cell_Stereo-seq/STDS0000058/stomics/E9.5_E1S1.MOSTA.h5ad")
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)

# %%
d.acquire_weight()
d.acquire_hotspot(cores=10)
d.acquire_density(cores=5)

# %%
d.find_clusters(n_svgs=500, n_svg_clusters=12)
d.svg_cluster.to_csv(f"results/E9.5_E1S1.MOSTA/svg_cluster.csv")
d.AI.sort_values(ascending=False).to_csv(f"results/E9.5_E1S1.MOSTA/AI.csv")

# %%
adjusted_coor_df = deepcopy(d.coordinate_df)
adjusted_coor_df["X"] -= (abs(min(d.coordinate_df["X"])) - 10)
adjusted_coor_df["Y"] += abs(min(d.coordinate_df["Y"]))
fig, ax = plt.subplots(figsize=(15, 10))
ax.axis("off")
ax.imshow(
    np.ones((
        int(max(adjusted_coor_df["X"])),
        int(max(adjusted_coor_df["Y"])),
        3,
    )))
sc = ax.scatter(
    adjusted_coor_df["X"],
    adjusted_coor_df["Y"],
    s=16,
    c=d.spot_type["type_1"],
    cmap="tab20",
)
legend = ax.legend(
    *sc.legend_elements(),
    bbox_to_anchor=(1, 1),
    title="Cluster",
    scatterpoints=3,
    ncol=2,
)
ax.add_artist(legend)
fig.savefig("results/E9.5_E1S1.MOSTA/typemap.1.jpg", bbox_inches="tight")

# %%
for gene in d.AI.sort_values(ascending=False)[:50].index:
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.imshow(
        np.ones((
            int(max(adjusted_coor_df["X"])),
            int(max(adjusted_coor_df["Y"])),
            3,
        )))
    ax.axis("off")
    sc = ax.scatter(
        adjusted_coor_df["X"],
        adjusted_coor_df["Y"],
        s=16,
        c=d.count_df[gene],
        cmap="Spectral_r",
    )
    ax.set_title(
        f"{gene}, in cluster {d.svg_cluster[gene]}. AI = {d.AI[gene]:.2f}")
    fig.colorbar(sc)
    fig.savefig(f"results/E9.5_E1S1.MOSTA/{gene}.jpg", bbox_inches="tight")
    plt.close(fig)