#! venv/bin/python
# -*- encoding: utf-8 -*-

# %%
from pathlib import Path

import src.utils.read_func as read_func
import src.utils.run_func as run_func
import src.utils.read_label as read_label

WORKDIR = read_func.WORKDIR

# %% mba
gsm = "GSM4459950"
sample = "18A"
count_df, coor_df = read_func.read_mba(gsm, sample)
result = run_func.run_spatialde(count_df, coor_df)
result.to_csv(
    Path.joinpath(WORKDIR, f"results/spatialde/{sample}-spatialde.csv"))

# %% als
sample = "GSM3399149_CN68_E1"
count_df, coor_df = read_func.read_als(sample)
result = run_func.run_spatialde(count_df, coor_df)
result.to_csv(
    Path.joinpath(WORKDIR, f"results/spatialde/{sample}-spatialde.csv"))

# %% stereo
sample = "E9.5_E1S1.MOSTA"
count_df, coor_df = read_func.read_stereo(sample)
result = run_func.run_spatialde(count_df, coor_df)
result.to_csv(
    Path.joinpath(WORKDIR, f"results/spatialde/{sample}-spatialde.csv"))

# %% stereo
sample = "Mouse_brain"
count_df, coor_df = read_func.read_stereo(sample)
if sample == "Mouse_brain":
    count_df = count_df[count_df.T.sum() > 1500]
    coor_df = coor_df.reindex(index=count_df.index)
result = run_func.run_spatialde(count_df, coor_df)
result.to_csv(
    Path.joinpath(WORKDIR, f"results/spatialde/{sample}-spatialde.csv"))

# %% DLPFC
sample = 151673
samples = [
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
]
for sample in samples:
    count_df, coor_df, _ = read_func.read_dlpfc(sample)
    coor_df.columns = ["X", "Y"]
    labels = read_label.read_dlpfc_label(sample)
    na_spots = labels[labels.isna()].index
    count_df = count_df.drop(index=na_spots)
    coor_df = coor_df.drop(index=na_spots)
    count_df = count_df.reindex(columns=count_df.columns[count_df.sum() > 0])
    result = run_func.run_spatialde(count_df, coor_df)
    result.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/spatialde/DLPFC-{sample}-spatialde.csv",
        ))
