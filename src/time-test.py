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
import time

import svgbit as sb

import src.utils.read_func as read_func
import src.utils.run_func as run_func

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

# %%
f = open("time1.txt", "a")
for idx in idx_full:
    d = sb.load_10X(f"data/mouse-brain-full/{idx}/outs")
    d = sb.filters.low_variance_filter(d)
    d = sb.filters.high_expression_filter(d)
    d = sb.normalizers.logcpm_normalizer(d)

    start_time = time.time()
    sb.run(d, cores=5)
    end_time = time.time()

    f.write(f"{idx},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
    f.write(f"{end_time - start_time}\n")

f.close()

# %%
f = open("time2.txt", "a")
for idx in idx_full:
    d = sb.load_10X(f"data/mouse-brain-full/{idx}/outs")
    d = sb.filters.low_variance_filter(d)
    d = sb.filters.high_expression_filter(d)

    start_time = time.time()
    run_func.run_somde(d.count_df, d.coordinate_df)
    end_time = time.time()

    f.write(f"{idx},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
    f.write(f"{end_time - start_time}\n")

f.close()

# %%
f = open("time3.txt", "a")
for idx in idx_full:
    d = sb.load_10X(f"data/mouse-brain-full/{idx}/outs")
    d = sb.filters.low_variance_filter(d)
    d = sb.filters.high_expression_filter(d)

    start_time = time.time()
    run_func.run_spatialde(d.count_df.sparse.to_dense(), d.coordinate_df)
    end_time = time.time()

    f.write(f"{idx},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
    f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
    f.write(f"{end_time - start_time}\n")

f.close()

# %%
idx = 151673
count_df, coor_df, *_ = read_func.read_dlpfc(idx)

f = open("time1.txt", "a")
d = sb.STDataset(count_df, coor_df)
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)
d = sb.normalizers.logcpm_normalizer(d)
start_time = time.time()
sb.run(d, cores=5)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()

f = open("time2.txt", "a")
start_time = time.time()
run_func.run_somde(count_df, coor_df)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()

f = open("time3.txt", "a")
start_time = time.time()
run_func.run_spatialde(count_df, coor_df)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()

# %%
idx = "Mouse_brain"
count_df, coor_df = read_func.read_stereo(idx)

# %%
f = open("time1.txt", "a")
d = sb.STDataset(count_df, coor_df)
d = sb.filters.low_variance_filter(d)
d = sb.filters.high_expression_filter(d)
d = sb.normalizers.logcpm_normalizer(d)
start_time = time.time()
sb.run(d, cores=5)
end_time = time.time()


start_time = time.time()
d.acquire_density(cores=1)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()

# %%
f = open("time2.txt", "a")
start_time = time.time()
run_func.run_somde(count_df, coor_df)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()

# %%
f = open("time3.txt", "a")
count_df_sde = count_df[count_df.T.sum() > 1500]
coor_df_sde = coor_df.reindex(index=count_df_sde.index)
start_time = time.time()
run_func.run_spatialde(count_df_sde, coor_df_sde)
end_time = time.time()
f.write(f"{idx},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(start_time))},")
f.write(f"{time.strftime('%m-%d %H:%M:%S', time.localtime(end_time))},")
f.write(f"{end_time - start_time}\n")
f.close()
