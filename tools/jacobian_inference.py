import numpy as np
import pandas as pd
import curve_interpolation

def traj_inference(adata, width=0.1, inc=0.05, nsim=10, frac=0.9): #see if strong dependen on any of the params
    # plain spliceJAC inference using the original datapoints without smoothing
    pst = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))
    ints = curve_interpolation.define_intervals(pst, width=width, inc=inc)

    tm, npoints = np.zeros(len(ints)), np.zeros(len(ints))
    max_eig, pos = np.zeros((len(ints), nsim)), np.zeros((len(ints), nsim))

    eig_list, jac_list = [], []

    # use the imputated unspliced and spliced counts
    U_data, S_data = adata.layers['Mu'], adata.layers['Ms']

    for i in range(len(ints)): #slice in pseudotime
        print('Running inference on pseudotime interval: ' + str(ints[i]))
        t1, t2 = ints[i][0], ints[i][1]
        sel_pst, sel_U, sel_S = pst[pst>t1], U_data[pst>t1], S_data[pst>t1]
        sel_U, sel_S = sel_U[sel_pst<t2], sel_S[sel_pst<t2]

        npoints[i] = sel_U.shape[0]
        tm[i] = (t1 + t2) / 2.

        int_mat = np.zeros((U_data.shape[1], U_data.shape[1]))

        # for each pseudotime point, run the inference multiple times (nsim)
        # each time, only a randomly selected fraction of the data is used
        for j in range(nsim):
            indices = np.sort(
                np.random.choice(np.arange(0, sel_U.shape[0], 1, dtype='int'), size=int(frac*sel_U.shape[0]), replace=False))
            B, C, G = parameter_regression(sel_U[indices], sel_S[indices])
            int_mat = int_mat + B

            J = construct_jac(B, G)
            w, v = np.linalg.eig(J)
            w = np.real(w)
            max_eig[i][j], pos[i][j] = np.amax(w), w[w > 0].size
        # average jacobian and spectrum
        int_mat = int_mat/nsim
        J = construct_jac(int_mat, G)
        w, v = np.linalg.eig(J)
        w = np.real(w)
        eig_list.append(w)
        jac_list.append(J)

        #instead of below need a different line defining the new name like the dictionary
    adata.uns['Jacobian'] = {'time': tm, 'largest': max_eig, 'positive': pos, 'spectrum':eig_list, 'jacobians':jac_list,
                             'npoints': npoints, 'pst_interval': ints,
                             'inference_params':{'width':width, 'inc':inc, 'nsim':nsim, 'frac':frac}}

def cell_capture(adata, width=0.1, inc=0.05, nsim=10, frac=0.9):
    # plain spliceJAC inference using the original datapoints without smoothing
    pst = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))
    ints = define_intervals(pst, width=width, inc=inc)

    tm, npoints = np.zeros(len(ints)), np.zeros(len(ints))
    cell_list, geneexp_list = [], []

    # use the imputated unspliced and spliced counts
    # U_data, S_data = adata.layers['Mu'], adata.layers['Ms']
    # Use the cell count
    Cell_data = adata.X

    for i in range(len(ints)):
        print('Running inference on pseudotime interval: ' + str(ints[i]))
        t1, t2 = ints[i][0], ints[i][1]

        if t1 < np.min(pst) or t2 > np.max(pst):
            # Print an error message or handle the case where t1 or t2 is out of range.
            print(f"Error: t1={t1} or t2={t2} is out of range of pseudotime values.")
            t2 = np.max(pst)
            print(f'Assigned new t2 = {t2}')
            C_window = Cell_data[(pst > t1) & (pst < t2), :]
            avg_Cellexp = np.mean(C_window, axis=0)
            geneexp_list.append(avg_Cellexp)
            print(f"Error fixed proceed as normal")
        else:
            C_window = Cell_data[(pst > t1) & (pst < t2), :]
            avg_Cellexp = np.mean(C_window, axis=0)
            geneexp_list.append(avg_Cellexp)

        # print(f'Cells in window:{C_window.shape[0]}\nMean expression: {avg_Cellexp}')

        # npoints[i] = C_window.shape[0]
        tm[i] = (t1 + t2) / 2.

        adata.uns['Cells_Captured'] = {'time': tm, 'gene_expression': geneexp_list, 'cells': cell_list,
                                       'pst_interval': ints,
                                       'inference_params': {'width': width, 'inc': inc, 'nsim': nsim, 'frac': frac}}
