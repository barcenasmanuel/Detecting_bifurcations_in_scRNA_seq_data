import numpy as np
from sklearn.linear_model import Ridge, LinearRegression, Lasso


def parameter_regression(U_data,
                         S_data,
                         method='Ridge',
                         alpha=1,
                         fit_int=True
                         ):
    """
        Run regression to infer spliced-unspliced interaction coefficients.

        Parameters:
            U_data (numpy.ndarray): Count matrix of unspliced counts.
            S_data (numpy.ndarray): Count matrix of spliced counts.
            method (str, optional): Regression method. Choose between 'Linear', 'Ridge', or 'Lasso'.
            alpha (float, optional): Regularization coefficient for Ridge and Lasso.
            fit_int (bool, optional): If True, set the fit_intercept parameter to True.

        Returns:
            mat (numpy.ndarray): Gene-gene interaction matrix.
            interc (numpy.ndarray): Intercept vector.
            degr (numpy.ndarray): Degradation coefficient vector.
    """
    assert method == 'Ridge' or method == 'Lasso' or method == 'Linear', "Please choose method='Ridge', 'Lasso' or 'Linear' "

    if method == 'Linear':
        reg = LinearRegression(fit_intercept=fit_int)
    elif method == 'Ridge':
        reg = Ridge(alpha=alpha, fit_intercept=fit_int)
    elif method == 'Lasso':
        reg = Lasso(alpha=alpha, fit_intercept=fit_int)

    ncell, ngene = U_data.shape
    mat = np.zeros((ngene, ngene))
    interc = np.zeros(ngene)
    degr = np.zeros(ngene)

    for i in range(ngene):
        S_use = np.delete(S_data, i, 1)
        reg.fit(S_use, U_data[:, i])
        coeffs = reg.coef_

        mat[i][0:i] = coeffs[0:i]
        mat[i][i + 1:] = coeffs[i:]
        interc[i] = reg.intercept_

        # fit spliced degradation rate - degradation rate sin spliceJAC are computed globally, so this step is optional
        reg_g = LinearRegression(fit_intercept=False)
        reg_g.fit(S_data[:, [i]], U_data[:, i])
        degr[i] = reg_g.coef_

    return mat, interc, degr


def estimate_degr(U_data, S_data
                  ):
    """
        Estimate degradation rates from unspliced and spliced count matrices.

        Parameters:
            U_data (numpy.ndarray): Count matrix of unspliced counts.
            S_data (numpy.ndarray): Count matrix of spliced counts.

        Returns:
            degr (numpy.ndarray): Degradation coefficient vector.
    """
    ncell, ngene = U_data.shape
    degr = np.zeros(ngene)

    # global fit using all cells since degradation rates are global
    for i in range(ngene):
        reg_g = LinearRegression(fit_intercept=False)
        reg_g.fit(S_data[:, [i]], U_data[:, i])
        degr[i] = reg_g.coef_
    return degr


def construct_jac(mat,
                  degr,
                  b=1
                  ):
    """
        Construct a Jacobian matrix given the gene-gene interactions and degradation rates.

        Parameters:
            mat (numpy.ndarray): Matrix of gene-gene interactions.
            degr (numpy.ndarray): Degradation coefficient vector.
            b (float, optional): Splicing rate constant. (Default value is 1)

        Returns:
            J (numpy.ndarray): Jacobian matrix.
    """
    ngene = mat.shape[0]

    jac1 = np.diag(- b *np.ones(ngene))   # unspliced-unspliced part
    jac2 = np.diag( b *np.ones(ngene))    # spliced-unspliced part
    jac3 = np.diag(-degr)               # spliced-spliced part

    J1 = np.concatenate([jac1, mat], axis=1)
    J2 = np.concatenate([jac2, jac3], axis=1)
    J = np.concatenate([J1, J2])

    return J
