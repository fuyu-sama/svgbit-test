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

import svgbit as sb

import src.utils.read_func as read_func
import src.utils.read_label as read_label
from src.utils.read_func import WORKDIR

# %%
samples = [
    151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673,
    151674, 151675, 151676
]
for sample in samples:
    count_df, coor_df, logcount_df = read_func.read_dlpfc(sample)
    labels = read_label.read_dlpfc_label(sample)
    na_spots = labels[labels.isna()].index
    d = sb.STDataset(logcount_df.drop(index=na_spots),
                     coor_df.drop(index=na_spots))
    d = sb.filters.low_variance_filter(d)

    sb.run(d, cores=10)

    write_dir = Path.joinpath(WORKDIR, f"results/DLPFC-{sample}")
    if not write_dir.exists():
        write_dir.mkdir()
    d.AI.sort_values(ascending=False).to_csv(Path.joinpath(
        write_dir, "AI.csv"))
    d.svg_cluster.to_csv(Path.joinpath(write_dir, "svg_cluster.csv"))
    sb.plot.spot_type_map(d, Path.joinpath(write_dir, "typemap.jpg"))
