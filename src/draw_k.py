#! venv/bin/python
# -*- encoding: utf-8 -*-

# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({"font.size": 18})


def perform_k_coor(sample):
    spearman_df = pd.read_csv(
        f"results/k_corr/{sample}-spearman.csv",
        index_col=0,
        header=0,
    )
    pearson_df = pd.read_csv(
        f"results/k_corr/{sample}-pearson.csv",
        index_col=0,
        header=0,
    )

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(
        spearman_df,
        ax=ax,
        annot=True,
        linewidth=1,
        linecolor="black",
        cmap="binary_r",
        vmax=-1,
        vmin=-2,
        cbar=False,
    )
    ax.set_title(sample)
    ax.set_xlabel("K")
    ax.set_ylabel("K")
    ax.invert_yaxis()
    fig.savefig(f"results/k_corr/{sample}-spearman.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(
        pearson_df,
        ax=ax,
        annot=True,
        linewidth=1,
        linecolor="black",
        cmap="binary_r",
        vmax=-1,
        vmin=-2,
        cbar=False,
    )
    ax.set_title(sample)
    ax.set_xlabel("K")
    ax.set_ylabel("K")
    ax.invert_yaxis()
    fig.savefig(f"results/k_corr/{sample}-pearson.svg")
    plt.close(fig)


# %%
sample = 151673
perform_k_coor(f"DLPFC-{sample}")

# %%
for sample in [
        "E135A", "E135B", "E155A", "E155B", "E165A", "E165B", "E175A1",
        "E175A2", "P0A1", "P0A2"
]:
    perform_k_coor(sample)
