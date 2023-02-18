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
gsm = "GSM4459950"
idx = "18A"
count_path = f"data/2020_SciAdv_MouseBrainAtlus/COUNT/{gsm}_expr_raw_counts_table_{idx}.tsv.gz"
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
    sep="\t",
    compression="gzip",
)

coor_path = f"data/2020_SciAdv_MouseBrainAtlus/COORDINATES/coordinates_{idx}.tsv"
coor_df = pd.read_csv(coor_path, index_col=0, header=0, sep="\t")

d = sb.STDataset(count_df, coor_df)
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d, 0.99)
d = sb.normalizers.logcpm_normalizer(d)

sb.run(d, cores=5, n_svg_clusters=10, n_svgs=500, k=8)
save_dir = Path(f"results/{idx}")
if not save_dir.exists():
    save_dir.mkdir()
d.svg_cluster.to_csv(f"results/{idx}/svg_cluster.csv")
d.AI.sort_values(ascending=False).to_csv(f"results/{idx}/AI.csv")
