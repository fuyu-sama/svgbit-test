#! venv/bin/python
# -*- encoding: utf-8 -*-

# %%
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

import svgbit as sb

plt.rcParams.update({"font.size": 18})

# %%
d = sb.load_anndata_h5(
    "data/2022_Cell_Stereo-seq/STDS0000058/stomics/Mouse_brain.h5ad")
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)

# %%
d.acquire_weight()
d.acquire_hotspot(cores=10)
d.acquire_density(cores=5)

# %%
d.find_clusters(n_svgs=500, n_svg_clusters=12)
d.AI.sort_values(ascending=False).to_csv("results/Mouse_brain/AI.csv")

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
fig.savefig("results/Mouse_brain/typemap.1.jpg", bbox_inches="tight")

# %%
# for gene in d.AI.sort_values(ascending=False)[:50].index:
for gene in ["Zbtb20"]:
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
    fig.savefig(f"results/Mouse_brain/{gene}.jpg", bbox_inches="tight")
    plt.close(fig)
