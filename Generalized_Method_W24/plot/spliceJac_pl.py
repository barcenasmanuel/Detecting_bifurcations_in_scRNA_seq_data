import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import scanpy as sc
import networkx as nx
import os

from ..tools import curve_interpolation

small_size = 6
medium_size = 8
large_size = 12


def high_expr_genes(adata, figname):
    fig = plt.figure(figsize=(4, 12))
    ax = plt.subplot(111)
    sc.pl.highest_expr_genes(adata, n_top=30, ax=ax, show=False)
    plt.tight_layout()
    plt.savefig(figname, format='pdf', dpi=300)


def pst_barplot(adata, groupby='time', order=None, figname=None):
    if order:
        keys = order
    else:
        keys = list(set(list(adata.obs[groupby])))

    pst_list = []
    for k in keys:
        pst_list.append(np.ndarray.flatten(np.asarray(adata[adata.obs['time'] == k].obs['dpt_pseudotime'])))

    plt.boxplot(pst_list)
    plt.xticks(np.arange(1, len(keys) + 1, 1), keys)
    plt.xlabel('Cell group')
    plt.ylabel('Pseudotime')

    if figname:
        plt.tight_layout()
        plt.savefig(figname, format='pdf', dpi=300)


def top_genes_pseudotime(adata, gene_list, figname, groupby='time', order=None):
    fig = plt.figure(figsize=(10.62, 10.62))
    ax_list = [331, 332, 333, 334, 335, 336, 337, 338, 339]
    for ax, gene in zip(ax_list, gene_list):
        one_gene(adata, gene, groupby=groupby, figname=None, order=order, ax=plt.subplot(ax))
    plt.tight_layout()
    plt.savefig(figname, format='pdf', dpi=300)


def one_gene(adata, gene, groupby='time', figname=None, order=None, ax=None):
    x = np.ndarray.flatten(np.asarray(adata.obs['dpt_pseudotime']))
    y_ms = adata[:, gene].layers['Ms']

    S_data = adata.uns['S_curves']
    pst, smoothed = S_data['avg_pst'], S_data[gene]

    states = np.asarray(adata.obs['time'])
    if order:
        keys = order
    else:
        keys = list(set(list(adata.obs[groupby])))

    if ax == None:
        fig = plt.figure(figsize=(4, 3))
        ax = plt.subplot(111)
    for k in keys:
        ax.scatter(x[states == k], y_ms[states == k], label=k, alpha=0.5)
    ax.plot(pst, smoothed, 'k-')
    ax.set_xlabel('Pseudotime')
    ax.set_ylabel(gene)
    ax.legend(fontsize=8)
    if figname:
        plt.tight_layout()
        plt.savefig(figname, format='pdf', dpi=300)


def plot_max_eig(adata, key='positive', figname=None, ax=None):
    assert key == 'positive' or key == 'max', "Choose between key=='positive' or key='max'"

    time, largest, positive = adata.uns['Jacobian']['time'], adata.uns['Jacobian']['largest'], \
        adata.uns['Jacobian']['positive']

    avg_pos, avg_max = np.mean(positive, axis=1), np.mean(largest, axis=1)
    std_pos, std_max = np.std(positive, axis=1), np.std(largest, axis=1)

    generate_fig = False
    if ax == None:
        fig = plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        generate_fig = True

    if key == 'max':
        ax.errorbar(time, avg_max, yerr=std_max, fmt='o')
        plt.plot([np.amin(time), np.amax(time)], [0, 0], 'r--')
        plt.ylabel('Largest eigenvalue')
    else:
        ax.errorbar(time, avg_pos, yerr=std_pos, fmt='o')
        plt.xlabel('Pseudotime')
        plt.ylabel('Positive eigenvalues')

    if generate_fig:
        plt.tight_layout()
        if figname == None:
            figname = 'traj_stability.pdf'
        plt.savefig(figname, format='pdf', dpi=300, bbox_inches='tight')


def plot_eigenspectrum(adata, point, title=True, figname=None, ax=None):
    time, eigen = adata.uns['Jacobian']['time'], adata.uns['Jacobian']['spectrum']

    generate_fig = False
    if ax == None:
        fig = plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        generate_fig = True
    ax.plot(np.flip(np.sort(eigen[point]))[0:10], 'o')
    plt.plot([-0.5, 9.5], [0, 0], 'r--')
    plt.xlim([-0.5, 9.5])
    plt.xlabel('Eigenvalue rank')
    plt.ylabel('Largest eigenvalues')
    if title:
        plt.title('pst=' + str(np.round(time[point], 2)))

    if generate_fig:
        plt.tight_layout()
        if figname == None:
            figname = 'spectrum_detail.pdf'
        plt.savefig(figname, format='pdf', dpi=300)


def plot_stability(adata, filename_p, rep=1, sizefig=3.54, tick_size=small_size, axis_size=medium_size, format='png', dpi=300):
    # Ensure that the target directory exists
    target_directory = f'figures/{filename_p}/'
    os.makedirs(target_directory, exist_ok=True)

    max_eig = adata.uns['Jacobian']['largest']
    tm = adata.uns['Jacobian']['time']
    pos = adata.uns['Jacobian']['positive']
    ints = adata.uns['Jacobian']['pst_interval']
    npoints = adata.uns['Jacobian']['npoints']

    avg_max, avg_pos = np.mean(max_eig, axis=1), np.mean(pos, axis=1)
    std_max, std_pos = np.std(max_eig, axis=1), np.std(pos, axis=1)

    smooth_tm, smooth_avg = curve_interpolation.smooth_curve(tm, avg_max, ints, rep=rep)
    smooth_tm, smooth_pos = curve_interpolation.smooth_curve(tm, avg_pos, ints, rep=rep)

    matplotlib.rc('xtick', labelsize=tick_size)
    matplotlib.rc('ytick', labelsize=tick_size)
    fig = plt.figure(figsize=(sizefig, 3*sizefig))  # figsize=(WIDTH_SIZE,HEIGHT_SIZE); figsize = (3.54,3.54)

    ax1 = plt.subplot(311)
    plt.errorbar(tm, avg_pos, yerr=std_pos, fmt='o')
    plt.plot(smooth_tm, smooth_pos, 'r--')
    plt.ylabel('Positive eigenvalues', fontsize=axis_size)
    plt.xlabel('Pseudotime', fontsize=axis_size)
    # plt.xticks([])

    ax2 = plt.subplot(312)
    plt.errorbar(tm, avg_max, yerr=std_max, fmt='o')
    plt.plot(smooth_tm, smooth_avg, 'r--')
    plt.ylabel('Largest eigenvalue', fontsize=axis_size)
    plt.xlabel('Pseudotime', fontsize=axis_size)
    # plt.xticks([])

    ax3 = plt.subplot(313)
    plt.plot(tm, npoints, 'o')
    plt.plot([0, 1], [50, 50], 'r--')
    plt.xlabel('Pseudotime', fontsize=axis_size)
    plt.ylabel('Number of observables', fontsize=axis_size)

    plt.tight_layout()

    plt.savefig(f'{target_directory}regression_{filename_p}.{format}', format=format, bbox_inches='tight', dpi=dpi)
