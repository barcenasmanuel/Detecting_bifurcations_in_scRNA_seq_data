import numpy as np
import networkx as nx
from sklearn.linear_model import Ridge, LinearRegression, Lasso


def parameter_regression(U_data,
                         S_data,
                         method='Ridge',
                         alpha=1,
                         fit_int=True
                         ):
    '''Run regression to infer spliced-unspliced interaction coefficients
    Parameters
    ----------
    U_data: `~numpy.ndarray`
        count matrix of unspliced counts
    S_data: `~numpy.ndarray`
        count matrix of spliced counts
    method: `str` (default: Ridge)
        regression method, choose between Linear, Ridge or Lasso
    alpha: `float` (default: 1)
        regularization coefficient for Ridge and Lasso
    fit_int: `Bool` (default: True)
        if True, set the fit_intercept parameter to True
    Returns
    -------
    mat: `~numpy.ndarray`
        gene-gene interaction matrix
    interc: `~numpy.ndarray`
        intercept vector
    degr: `~numpy.ndarray`
        degradation coefficient vector
    '''
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
    '''Construct a Jacobian matrix given the gene-gene interactions and degradation rates
    Parameters
    ----------
    mat: `~numpy.ndarray`
        matrix of gene-gene interactions
    degr: `~numpy.ndarray`
        degradation coefficient vector
    b: `float` (default: 1)
        splicing rate constant
    Returns
    -------
    J: Jacobian matrix
    '''
    ngene = mat.shape[0]

    jac1 = np.diag(- b *np.ones(ngene))   # unspliced-unspliced part
    jac2 = np.diag( b *np.ones(ngene))    # spliced-unspliced part
    jac3 = np.diag(-degr)               # spliced-spliced part

    J1 = np.concatenate([jac1, mat], axis=1)
    J2 = np.concatenate([jac2, jac3], axis=1)
    J = np.concatenate([J1, J2])

    return J


def calc_grn(adata, genes=None, index=0, weight_quantile=.5):
    if genes == None:
        genes = list(adata.var_names)

    n = len(genes)
    A = adata.uns['Jacobian']['jacobians'][index][0:n, n:].copy().T

    q_pos = np.quantile(A[A > 0], weight_quantile)
    q_neg = np.quantile(A[A < 0], 1 - weight_quantile)
    A[(A > q_neg) & (A < q_pos)] = 0

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph)
    nx.relabel_nodes(G, dict(zip(range(len(genes)), genes)), copy=False)

    return G


def grn_geneweight_exp(adata, genes=None, index=0, weight_quantile=.5):
    if genes is None:
        # print(f'No genes detect...using adata.var_names')
        genes = list(adata.var_names)
        # print (f'Genes (lenght= {len(genes)}): {genes}')

    n = len(genes)
    A = adata.uns['Jacobian']['jacobians'][index][0:n, n:].copy().T

    # Calculate the average gene expression for each edge
    gene_expression_2d = adata.uns['Cells_Captured']['gene_expression'][index]
    gene_expression = gene_expression_2d.flatten()

    # print (f'Gene_expression {index}: {gene_expression}') #see difference

    q_pos = np.quantile(A[A > 0], weight_quantile)
    q_neg = np.quantile(A[A < 0], 1 - weight_quantile)
    A[(A > q_neg) & (A < q_pos)] = 0

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph)
    nx.relabel_nodes(G, dict(zip(range(len(genes)), genes)), copy=False)

    for i, gene1 in enumerate(genes):
        for j, gene2 in enumerate(genes):
            if i != j:
                if G.has_edge(gene1, gene2):
                    default_weight = G[gene1][gene2].get('weight', 1.0)
                    weight = default_weight * gene_expression[i]
                    G[gene1][gene2]['weight'] = weight
                else:
                    # Handle it?
                    pass

    return G

def G_listgen(adata, jac_list):
    G_list = []
    print(f'Running comm_tools.G_listgen')
    for i, jac_value in enumerate(jac_list):
        G = calc_grn(adata, index=i,weight_quantile=0.9) #wq=0.9

        G_list.append(G)
    #print('Done')
    return G_list


def G_listgen_geneexp(adata, jac_list, geneexp_list):
    print(f'Running comm_tools.G_listgen_geneexp')
    G_list = []
    for i, jac_value in enumerate(jac_list):
        G = grn_geneweight_expression2(adata, index=i, weight_quantile=0.9)  # wq=0.9
        G_list.append(G)

    return G_list
