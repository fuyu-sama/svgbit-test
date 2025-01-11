#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
from pathlib import Path

import matplotlib.pyplot as plt

import src.utils.read_func as read_func
from src.utils.read_func import WORKDIR

# %%
sample = "E9.5_E1S1.MOSTA"
reads = read_func.read_stereo(sample)
genes = []
with open("genes.txt") as f:
    for line in f:
        line = line.strip()
        genes.append(line)
write_dir = Path.joinpath(WORKDIR, "gene-expression", sample)
if not write_dir.exists():
    write_dir.mkdir()

# %%
for gene in genes:
    fig, ax = plt.subplots(figsize=(10, 10))
    scp = ax.scatter(
        reads.coor_df["X"],
        reads.coor_df["Y"],
        s=4,
        c=reads.count_df[gene],
        cmap="copper_r",
        alpha=0.7,
        vmax=reads.count_df[gene].quantile(0.99) + 1,
        vmin=-(reads.count_df[gene].quantile(0.99) + 1) / 4,
    )
    ax.axis("off")
    fig.colorbar(scp, ax=ax)
    fig.savefig(Path.joinpath(write_dir, f"{gene}.jpg"))
    plt.close(fig)
