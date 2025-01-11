#! venv/bin/python
# -*- encoding: utf-8 -*-

# %%
import pickle
from copy import deepcopy
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import svgbit as sb

import src.utils.read_func as read_func
from src.utils.read_func import WORKDIR

# %% dlpfc
samples = [
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
]
for sample in samples:
    reads = read_func.read_dlpfc(sample)
    d = sb.STDataset(reads.logcount_df, reads.coor_df)
    d = sb.filters.low_variance_filter(d)

    sb.run(d, cores=10)

    write_dir = Path.joinpath(WORKDIR, f"results/svgbit/DLPFC-{sample}")
    if not write_dir.exists():
        write_dir.mkdir()
    d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
        write_dir, "AI.csv"))
    d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
    sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
    with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
        pickle.dump(d, f)

# %% GBM
samples = ["21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2"]
for sample in samples:
    reads = read_func.read_gbm(sample)
    d = sb.STDataset(reads.count_df, reads.coor_df)
    d = sb.filters.low_variance_filter(d)
    d = sb.normalizers.logcpm_normalizer(d)

    sb.run(d, cores=10)

    write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
    if not write_dir.exists():
        write_dir.mkdir()
    d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
        write_dir, "AI.csv"))
    d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
    sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
    with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
        pickle.dump(d, f)

# %%
sample = "E9.5_E1S1.MOSTA"
write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
if not write_dir.exists():
    write_dir.mkdir()

d = sb.load_anndata_h5(
    f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)

d.acquire_weight()
d.acquire_hotspot(cores=10)
d.acquire_density(cores=5)
d.find_clusters(n_svgs=500, n_svg_clusters=12)

d.AI.sort_values(ascending=False).to_csv(Path.joinpath(write_dir, "AI.csv"))
with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
    pickle.dump(d, f)

# %% stereo
sample = "Mouse_brain"
write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
if not write_dir.exists():
    write_dir.mkdir()

d = sb.load_anndata_h5(
    f"data/2022_Cell_Stereo-seq/STDS0000058/stomics/{sample}.h5ad")
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)

d.acquire_weight()
d.acquire_hotspot(cores=10)
d.acquire_density(cores=2)
d.find_clusters(n_svgs=500, n_svg_clusters=12)

d.AI.sort_values(ascending=False).to_csv(Path.joinpath(write_dir, "AI.csv"))
with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
    pickle.dump(d, f)

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
fig.savefig(Path.joinpath(write_dir, "typemap.jpg"), bbox_inches="tight")

# %% iimpact
samples = [
    "human_breast_cancer", "human_ovarian_cancer", "human_prostate_cancer"
]
for sample in samples:
    reads = read_func.read_iimpact(sample)
    d = sb.STDataset(reads.count_df, reads.coor_df)
    d = sb.filters.low_variance_filter(d)
    d = sb.normalizers.logcpm_normalizer(d)

    sb.run(d, cores=10)

    write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
    if not write_dir.exists():
        write_dir.mkdir()
    d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
        write_dir, "AI.csv"))
    d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
    sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
    with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
        pickle.dump(d, f)

# %% mba
for sample in read_func.mba_files.keys():
    reads = read_func.read_mba(sample)
    try:
        d = sb.STDataset(reads.count_df, reads.coor_df)
        d = sb.filters.low_variance_filter(d)
        d = sb.normalizers.logcpm_normalizer(d)

        sb.run(d, cores=10)

        write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
        if not write_dir.exists():
            write_dir.mkdir()
        d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
            write_dir, "AI.csv"))
        d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
        sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
        with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
            pickle.dump(d, f)
    except:
        continue

# %%
samples = [
    "E135A", "E135B", "E155A", "E155B", "E165A", "E165B", "E175A1", "E175A2",
    "E175B", "P0A1", "P0A2", "P0B"
]
for sample in samples:
    reads = read_func.read_self(sample)
    d = sb.STDataset(reads.count_df, reads.coor_df)
    d = sb.filters.low_variance_filter(d)
    d = sb.normalizers.logcpm_normalizer(d)

    sb.run(d, cores=10)

    write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{sample}")
    if not write_dir.exists():
        write_dir.mkdir()
    d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
        write_dir, "AI.csv"))
    d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
    sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
    with open(Path.joinpath(write_dir, f"{sample}.pkl"), "wb") as f:
        pickle.dump(d, f)
