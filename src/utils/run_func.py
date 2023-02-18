import numpy as np

import NaiveDE
import SpatialDE

from somde import SomNode

import svgbit as sb


def run_spatialde(count_df, coor_df):
    coor_df = coor_df.copy()
    coor_df["total_counts"] = count_df.T.sum()
    norm_expr = NaiveDE.stabilize(count_df.T).T
    resid_expr = NaiveDE.regress_out(
        coor_df,
        norm_expr.T,
        'np.log(total_counts)',
    ).T
    results = SpatialDE.run(coor_df[["X", "Y"]], resid_expr).sort_values(by="qval")

    return results


def run_somde(count_df, coor_df):
    # count_df: spots * genes
    coor_df = coor_df.copy()
    coor_df["total_count"] = count_df.T.sum()
    som = SomNode(coor_df[["X", "Y"]].values.astype(np.float32), 20)
    ndf, ninfo = som.mtx(count_df.T)
    nres = som.norm()
    result, SVnum = som.run()
    result = result.sort_values(by="qval")

    return result


def run_svgbit(count_df, coor_df):
    d = sb.STDataset(count_df, coor_df)
    sb.run(d)

    return d.AI.sort_values(ascending=False)
