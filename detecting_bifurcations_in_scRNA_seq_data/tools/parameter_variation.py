import numpy as np
from typing import List
from ..tools import curve_interpolation
from ..tools import grn_inference


def par_traj_inference(adata, width=0.1, inc=0.05, nsim=10, frac=0.9, method='Ridge', alpha=1, fit_int=True, b=1):
    # plain spliceJAC inference using the original datapoints without smoothing
    pst = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))
    ints = curve_interpolation.define_intervals(pst, width=width, inc=inc)

    tm, npoints = np.zeros(len(ints)), np.zeros(len(ints))
    max_eig, pos = np.zeros((len(ints), nsim)), np.zeros((len(ints), nsim))

    eig_list, jac_list = [], []

    # use the imputated unspliced and spliced counts
    U_data, S_data = adata.layers['Mu'], adata.layers['Ms']

    for i in range(len(ints)):
        print('Running inference on pseudotime interval: ' + str(ints[i]))
        t1, t2 = ints[i][0], ints[i][1]
        sel_pst, sel_U, sel_S = pst[pst > t1], U_data[pst > t1], S_data[pst > t1]
        sel_U, sel_S = sel_U[sel_pst < t2], sel_S[sel_pst < t2]

        npoints[i] = sel_U.shape[0]
        tm[i] = (t1 + t2) / 2.

        int_mat = np.zeros((U_data.shape[1], U_data.shape[1]))

        # for each pseudotime point, run the inference multiple times (nsim)
        # each time, only a randomly selected fraction of the data is used
        for j in range(nsim):
            indices = np.sort(
                np.random.choice(np.arange(0, sel_U.shape[0], 1, dtype='int'), size=int(frac * sel_U.shape[0]),
                                 replace=False))
            B, C, G = grn_inference.parameter_regression(sel_U[indices], sel_S[indices], method=method,
                                                               alpha=alpha, fit_int=fit_int)
            int_mat = int_mat + B

            J = grn_inference.construct_jac(B, G, b=b)
            w, v = np.linalg.eig(J)
            w = np.real(w)
            max_eig[i][j], pos[i][j] = np.amax(w), w[w > 0].size
        # average jacobian and spectrum
        int_mat = int_mat / nsim
        J = grn_inference.construct_jac(int_mat, G, b=b)
        w, v = np.linalg.eig(J)
        w = np.real(w)
        eig_list.append(w)
        jac_list.append(J)

    result_dict = {
        'time': tm,
        'largest': max_eig,
        'positive': pos,
        'spectrum': eig_list,
        'jacobians': jac_list,
        'npoints': npoints,
        'pst_interval': ints,
        'inference_params': {'width': width, 'inc': inc, 'nsim': nsim, 'frac': frac}
    }

    return result_dict


def par_width_var(processed_adata, wsweep: List[float] = None, method='Ridge', alpha=1, fit_int=True, b=1):
    '''
    :param processed_adata:
    :param wsweep:
    :param method:
    :param alpha:
    :param fit_int:
    :param b:
    :return:
    '''
    processed_adata.uns['var_par'] = {}
    # Default wsweep value if not provided
    default_wsweep = np.linspace(0.05, 0.2, num=4)
    default_wsweep = np.around(default_wsweep, 4)

    # Use provided wsweep or default
    wsweep=wsweep if wsweep is not None else default_wsweep

    # Check if wsweep is a vector with exactly 4 values
    if len(wsweep) != 4:
        raise ValueError("wsweep must be a vector with exactly 4 values.")

    for value in wsweep:
        result_dict = par_traj_inference(processed_adata, width=value, method=method, alpha=alpha, fit_int=fit_int, b=b)
        processed_adata.uns['var_par'][f'width={value}'] = result_dict


def par_inc_var(processed_adata, isweep: List[float] = None, method='Ridge', alpha=1, fit_int=True, b=1):
    processed_adata.uns['var_par'] = {}
    # Default isweep values if not provided
    default_isweep = np.linspace(0.005, 0.05, num=4)
    default_isweep = np.around(default_isweep, 4)

    # Use provided isweep or default
    isweep = isweep if isweep is not None else default_isweep

    # Check if isweep is a vector with exactly 4 values
    if len(isweep) != 4:
        raise ValueError("isweep must be a vector with exactly 4 values.")

    for value in isweep:
        result_dict = par_traj_inference(processed_adata, inc=value, method=method, alpha=alpha, fit_int=fit_int, b=b)
        processed_adata.uns['var_par'][f'inc={value}'] = result_dict


def par_nsim_var(processed_adata, nsweep: List[float] = None, method='Ridge', alpha=1, fit_int=True, b=1):
    processed_adata.uns['var_par'] = {}
    # Default nsweep value if not provided
    default_nsweep = np.linspace(10, 40, num=4)
    default_nsweep = np.around(default_nsweep, 0)

    # Use provided nsweep or default
    nsweep = nsweep if nsweep is not None else default_nsweep

    # Check if nsweep is a vector with exactly 4 values
    if len(nsweep) != 4:
        raise ValueError("nsweep must be a vector with exactly 4 values.")

    for value in nsweep:
        result_dict = par_traj_inference(processed_adata, nsim=value, method=method, alpha=alpha, fit_int=fit_int, b=b)
        processed_adata.uns['var_par'][f'nsim={value}'] = result_dict


def par_frac_var(processed_adata, fsweep: List[float] = None, method='Ridge', alpha=1, fit_int=True, b=1):
    processed_adata.uns['var_par'] = {}
    # Default fsweep value if not provided
    default_fsweep = np.linspace(0.3, 0.9, num=4)
    default_fsweep = np.around(default_fsweep, 4)

    # Use provided fsweep or default
    fsweep = fsweep if fsweep is not None else default_fsweep

    # Check if fsweep is a vector with exactly 4 values
    if len(fsweep) != 4:
        raise ValueError("fsweep must be a vector with exactly 4 values.")

    for value in fsweep:
        result_dict = par_traj_inference(processed_adata, frac=value, method=method, alpha=alpha, fit_int=fit_int, b=b)
        processed_adata.uns['var_par'][f'frac={value}'] = result_dict
