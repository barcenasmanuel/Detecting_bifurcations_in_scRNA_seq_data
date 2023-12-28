'''
Tools for curve interpolation
'''
import numpy as np
import pandas as pd


def define_intervals(pst, width=0.1, inc=0.001):
    start, end = np.amin(pst), np.amax(pst)
    w = width*(end-start)
    endpoint = start+w
    ints = []
    ints.append([start, start+w])
    while(endpoint<end):
        i1, i2 = ints[-1] if ints else (start, start + w)
        ints.append([i1+inc, i2+inc])
        endpoint = ints[-1][1]
    return ints


def smooth_curve(pst, counts, intervals, rep=1):
    ref_pst, ref_counts = pst, counts
    for j in range(rep):
        avg_pst, avg_count = np.zeros(len(intervals)), np.zeros(len(intervals))

        for i in range(len(intervals)):
            t1, t2 = intervals[i][0], intervals[i][1]
            sel_pst, sel_counts = ref_pst[ref_pst > t1], ref_counts[ref_pst > t1]
            sel_pst, sel_counts = sel_pst[sel_pst < t2], sel_counts[sel_pst < t2]
            avg_pst[i], avg_count[i] = (t1 + t2) / 2., np.mean(sel_counts)
        ref_pst, ref_counts = avg_pst, avg_count

    return avg_pst, avg_count


def curve_average(y, rep=1):
    ref_y = y
    for j in range(rep):
        avg_y = y.copy()

        for i in range(1, ref_y.size-1):
            avg_y[i] = (ref_y[i-1]+ref_y[i+1])/2.
        ref_y = avg_y

    return avg_y


def smooth_counts(adata, width=0.1, inc=0.001, rep=1):
    geneset = list(adata.var_names)
    pst = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))

    ints = define_intervals(pst, width=width, inc=inc)

    U_dict, S_dict = {}, {}

    for g in geneset:
        dat_U, dat_S = adata[:, g].layers['Mu'], adata[:, g].layers['Ms']
        # unspliced
        counts = np.ndarray.flatten(dat_U)
        avg_pst, avg_count = smooth_curve(pst, counts, ints)
        smooth_count = curve_average(avg_count, rep=rep)
        U_dict[g] = smooth_count
        # spliced
        counts = np.ndarray.flatten(dat_S)
        avg_pst, avg_count = smooth_curve(pst, counts, ints)
        smooth_count = curve_average(avg_count, rep=rep)
        S_dict[g] = smooth_count

    U_dict['avg_pst'], S_dict['avg_pst'] = avg_pst, avg_pst
    U, S = pd.DataFrame.from_dict(U_dict), pd.DataFrame.from_dict(S_dict)

    adata.uns['U_curves'] = U
    adata.uns['S_curves'] = S
