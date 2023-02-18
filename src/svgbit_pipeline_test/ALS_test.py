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
from PIL import Image

import pandas as pd
import matplotlib.pyplot as plt

import svgbit as sb

Image.MAX_IMAGE_PIXELS = None
plt.rcParams.update({"font.size": 18})

# %%
sample = "GSM3399149_CN68_E1"
sample_short = "_".join(sample.split("_")[1:])
count_path = f"data/2019_Science_ALS/COUNT/{sample}_stdata_aligned_counts_IDs.txt.gz"
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
    sep="\t",
    compression="gzip",
).T

x_coor = []
y_coor = []
for i in count_df.index:
    x_coor.append(float(i.split("_")[0]))
    y_coor.append(float(i.split("_")[1]))
coor_df = pd.DataFrame((x_coor, y_coor)).T
coor_df.index = count_df.index
coor_df.columns = [["X", "Y"]]

d = sb.STDataset(count_df, coor_df)
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d, 0.99)
d = sb.normalizers.logcpm_normalizer(d)

sb.run(d, cores=10, n_svg_clusters=10, n_svgs=500)
save_dir = Path(f"results/{sample}")
if not save_dir.exists():
    save_dir.mkdir()
d.svg_cluster.to_csv(f"results/{sample}/svg_cluster.csv")
d.AI.sort_values(ascending=False).to_csv(f"results/{sample}/AI.csv")
