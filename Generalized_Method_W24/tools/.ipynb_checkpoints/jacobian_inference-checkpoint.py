import numpy as np

from ..tools import curve_interpolation
from ..tools import grn_inference


def traj_inference(adata, pst_key='dpt_pseudotime', method='Ridge', alpha=1, fit_int=True, width=0.1, inc=0.05, nsim=10, frac=0.9, b=1):
    """
        Perform plain spliceJAC inference using the original datapoints without smoothing.

        Parameters:
            adata (AnnData): Annotated data object.
            pst_key (str): Key in `adata.obs` representing pseudotime.
            method (str, optional): Regression method. Choose between 'Linear', 'Ridge', or 'Lasso'.
            alpha (float, optional): Regularization coefficient for Ridge and Lasso.
            fit_int (bool, optional): If True, set the fit_intercept parameter to True.
            width (float, optional): Width of pseudotime intervals for inference.
            inc (float, optional): Increment for defining pseudotime intervals.
            nsim (int, optional): Number of simulations for each pseudotime interval.
            frac (float, optional): Fraction of data used for each simulation.
            b (float, optional): Splicing rate constant.

        Returns:
            None
    """
    pst = np.ndarray.flatten(np.asarray(adata.obs[pst_key]))
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
                np.random.choice(np.arange(0, sel_U.shape[0], 1, dtype='int'), size=int(frac*sel_U.shape[0]),
                                 replace=False))
            B, C, G = grn_inference.parameter_regression(sel_U[indices], sel_S[indices], method=method,
                                                               alpha=alpha, fit_int=fit_int)
            int_mat = int_mat + B

            J = grn_inference.construct_jac(B, G, b=b)
            w, v = np.linalg.eig(J)
            w = np.real(w)
            max_eig[i][j], pos[i][j] = np.amax(w), w[w > 0].size
        # average jacobian and spectrum
        int_mat = int_mat/nsim
        J = grn_inference.construct_jac(int_mat, G, b=b)
        w, v = np.linalg.eig(J)
        w = np.real(w)
        eig_list.append(w)
        jac_list.append(J)

    adata.uns['Jacobian'] = {'time': tm, 'largest': max_eig, 'positive': pos, 'spectrum':eig_list, 'jacobians':jac_list,
                             'npoints': npoints, 'pst_interval': ints,
                             'inference_params':{'width':width, 'inc':inc, 'nsim':nsim, 'frac':frac}}


def cell_capture(adata, filename, width=0.1, inc=0.05, nsim=10, frac=0.9):
    """
        Capture cells in the window of plain spliceJAC inference using the original datapoints without smoothing.

        Parameters:
            adata (AnnData): Annotated data object.
            filename (str): Base filename for saving the resulting data.
            width (float, optional): Width of pseudotime intervals for cell capture.
            inc (float, optional): Increment for defining pseudotime intervals.
            nsim (int, optional): Number of simulations for each pseudotime interval.
            frac (float, optional): Fraction of data used for each simulation.

        Returns:
            None
    """
    pst = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))
    ints = curve_interpolation.define_intervals(pst, width=width, inc=inc)

    tm, npoints = np.zeros(len(ints)), np.zeros(len(ints))
    cell_list, geneexp_list = [], []

    # use the imputated unspliced and spliced counts
    Cell_data = adata.X

    for i in range(len(ints)):
        print('Running inference for cell capture on pseudotime interval: ' + str(ints[i]))
        t1, t2 = ints[i][0], ints[i][1]

        if t1 < np.min(pst) or t2 > np.max(pst):
            # Print an error message or handle the case where t1 or t2 is out of range.
            print(f"Error: t1={t1} or t2={t2} is out of range of pseudotime values.")
            t2 = np.max(pst)
            print(f'Assigned new t2 = {t2}')
            C_window = Cell_data[(pst > t1) & (pst < t2), :]
            avg_Cellexp = np.mean(C_window, axis=0)
            geneexp_list.append(avg_Cellexp)
            C_index = np.ndarray.flatten(np.argwhere((pst > t1) & (pst < t2)))
            cell_list.append(C_index)
            print(f"Error fixed proceed as normal")
        else:
            C_window = Cell_data[(pst > t1) & (pst < t2), :]
            avg_Cellexp = np.mean(C_window, axis=0)
            geneexp_list.append(avg_Cellexp)
            C_index = np.ndarray.flatten(np.argwhere((pst > t1) & (pst < t2)))
            cell_list.append(C_index)

        tm[i] = (t1 + t2) / 2.

    np.save(f'cell_indices_{filename}.npy', cell_list)
    adata.uns['Cells_Captured'] = {'time': tm, 'gene_expression': geneexp_list, 'pst_interval': ints,
                                   'inference_params': {'width': width, 'inc': inc, 'nsim': nsim, 'frac': frac}}


def create_weights_geneexpress(processed_adata, filename, method='Ridge', alpha=1, fit_int=True,
                               width=0.1, inc=0.05, nsim=10, frac=0.9, b=1):
    """
        Create weights for gene expression using plain spliceJAC inference and cell capture.

        Parameters:
            processed_adata (AnnData): Processed annotated data object.
            filename (str): Base filename for saving the resulting data.
            method (str, optional): Regression method. Choose between 'Linear', 'Ridge', or 'Lasso'.
            alpha (float, optional): Regularization coefficient for Ridge and Lasso.
            fit_int (bool, optional): If True, set the fit_intercept parameter to True.
            width (float, optional): Width of pseudotime intervals for inference and cell capture.
            inc (float, optional): Increment for defining pseudotime intervals.
            nsim (int, optional): Number of simulations for each pseudotime interval.
            frac (float, optional): Fraction of data used for each simulation.
            b (float, optional): Splicing rate constant.

        Returns:
            None
    """
    traj_inference(processed_adata, method=method, alpha=alpha, fit_int=fit_int, width=width,
                   inc=inc, nsim=nsim, frac=frac, b=b)
    cell_capture(processed_adata,filename, width=width, inc=inc, nsim=nsim, frac=frac)
    processed_adata.write_h5ad('gene_exp_'+filename+'.h5ad')
