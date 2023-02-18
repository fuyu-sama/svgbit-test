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

import matplotlib.pyplot as plt

import svgbit as sb

Image.MAX_IMAGE_PIXELS = None
plt.rcParams.update({"font.size": 18})

# %%
sample = "V1_Breast_Cancer_Block_A_Section_1"

d = sb.load_10X(f"data/10XGenomics/{sample}/outs")
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d, 0.99)
d = sb.normalizers.logcpm_normalizer(d)
sb.run(d, cores=6, n_svg_clusters=10, n_svgs=500)
save_dir = Path(f"results/{sample}")
he_path = Path(f"data/10XGenomics/{sample}/HE/{sample}_image.tif")
if not he_path.exists():
    he_path = Path("data/10XGenomics/{sample}/HE/{sample}_image.jpg")
if not save_dir.exists():
    save_dir.mkdir()
d.svg_cluster.to_csv(f"results/{sample}/svg_cluster.csv")
d.AI.sort_values(ascending=False).to_csv(f"results/{sample}/AI.csv")

# %%
he_image = Image.open(he_path)
fig, ax = plt.subplots(figsize=(10, 10))
ax.axis("off")
ax.imshow(he_image)
sc = ax.scatter(
    d.coordinate_df["X"],
    d.coordinate_df["Y"],
    s=16,
    c=d.spot_type["type_1"],
    cmap="tab20",
)
legend = ax.legend(
    *sc.legend_elements(),
    bbox_to_anchor=(1.2, 1),
    title="Cluster",
    scatterpoints=3,
)
ax.add_artist(legend)
plt.tight_layout()
fig.savefig(Path.joinpath(save_dir, "typemap.1.jpg"), bbox_inches="tight")

# %%
for gene in d.AI.sort_values(ascending=False)[:50].index:
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("off")
    ax.imshow(he_image)
    sc = ax.scatter(
        d.coordinate_df["X"],
        d.coordinate_df["Y"],
        s=16,
        c=d.count_df[gene],
        cmap="Spectral_r",
    )
    ax.set_title(
        f"{gene}, in cluster {d.svg_cluster[gene]}. AI = {d.AI[gene]:.2f}")
    fig.colorbar(sc)
    fig.savefig(Path.joinpath(save_dir, f"{gene}.jpg"), bbox_inches="tight")
    plt.close(fig)
