#! venv/bin/python
# -*- encoding: utf-8 -*-

# %%
from pathlib import Path

import src.utils.read_func as read_func
import src.utils.run_func as run_func
from src.utils.read_func import WORKDIR

write_dir = Path.joinpath(WORKDIR, "results", "SOMDE")
if not write_dir.exists():
    write_dir.mkdir()

# %% DLPFC
samples = [
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
]
for sample in samples:
    reads = read_func.read_dlpfc(sample)
    result = run_func.run_somde(reads.count_df, reads.coor_df)
    result.to_csv(Path.joinpath(write_dir, f"DLPFC-{sample}.csv"))

# %% stereo
sample = "E9.5_E1S1.MOSTA"
reads = read_func.read_stereo(sample)
result = run_func.run_somde(reads.count_df, reads.coor_df)
result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))

# %% stereo
sample = "Mouse_brain"
reads = read_func.read_stereo(sample)
result = run_func.run_somde(reads.count_df, reads.coor_df)
result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))

# %% iimpact
samples = [
    "human_breast_cancer", "human_ovarian_cancer", "human_prostate_cancer"
]
for sample in samples:
    reads = read_func.read_iimpact(sample)
    result = run_func.run_somde(reads.count_df, reads.coor_df)
    result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))

# %% mba
for sample in read_func.mba_files.keys():
    reads = read_func.read_mba(sample)
    try:
        result = run_func.run_somde(reads.count_df, reads.coor_df)
        result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))
    except:
        continue

# %% GBM
samples = ["21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2"]
for sample in samples:
    reads = read_func.read_gbm(sample)
    result = run_func.run_somde(
        reads.count_df.sparse.to_dense(),
        reads.coor_df,
    )
    result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))

# %%
samples = [
    "E135A", "E135B", "E155A", "E155B", "E165A", "E165B", "E175A1", "E175A2",
    "E175B", "P0A1", "P0A2", "P0B"
]
for sample in samples:
    reads = read_func.read_self(sample)
    result = run_func.run_somde(
        reads.count_df.sparse.to_dense(),
        reads.coor_df,
    )
    result.to_csv(Path.joinpath(write_dir, f"{sample}.csv"))
