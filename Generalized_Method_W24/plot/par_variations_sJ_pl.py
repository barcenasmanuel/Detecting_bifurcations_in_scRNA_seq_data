import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import scanpy as sc
import os
from typing import List

from ..tools import curve_interpolation

small_size = 6
medium_size = 8
large_size = 12


def plot_traj_widthvar(processed_adata, filename=None, wsweep: List[float] = None, fig_width=5, fig_height=7,
                       legend_size=medium_size, file_format='pdf'):
    # Ensure that the target directory exists
    if filename is not None:
        target_directory = f'figures/{filename}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Warning: filename is None. Saving to 'figures/' directory.")

    # Default wsweep value if not provided
    default_wsweep = np.linspace(0.05, 0.2, num=4)
    default_wsweep = np.around(default_wsweep, 4)

    # Use provided wsweep or default
    wsweep = wsweep if wsweep is not None else default_wsweep

    # Check if wsweep is a vector with exactly 4 values
    if len(wsweep) != 4:
        raise ValueError("wsweep must be a vector with exactly 4 values.")

    # create dictionaries to store values for plotting
    max_eig_dc = {}
    tm_dc = {}
    pos_dc = {}
    ints_dc = {}
    npoints_dc = {}
    avg_max_dc = {}
    avg_pos_dc = {}
    std_max_dc = {}
    std_pos_dc = {}

    smooth_tm_dc = {}
    smooth_avg_dc = {}
    smooth_pos_dc = {}

    for i, value in enumerate(wsweep):
        key1 = "max_eig" + str(i + 1)
        key2 = "tm" + str(i + 1)
        key3 = "pos" + str(i + 1)
        key4 = "ints" + str(i + 1)
        key5 = "npoints" + str(i + 1)
        key6 = "avg_max" + str(i + 1)
        key7 = "avg_pos" + str(i + 1)
        key8 = "std_max" + str(i + 1)
        key9 = "std_pos" + str(i + 1)
        key10 = "smooth_tm" + str(i + 1)
        key11 = "smooth_avg" + str(i + 1)
        key12 = "smooth_pos" + str(i + 1)

        max_eig_dc[key1] = processed_adata.uns['var_par']['width=%s' % value]['largest']
        tm_dc[key2] = processed_adata.uns['var_par']['width=%s' % value]['time']
        pos_dc[key3] = processed_adata.uns['var_par']['width=%s' % value]['positive']
        ints_dc[key4] = processed_adata.uns['var_par']['width=%s' % value]['pst_interval']
        npoints_dc[key5] = processed_adata.uns['var_par']['width=%s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = curve_interpolation.smooth_curve(tm_dc[key2], avg_max_dc[key6],
                                                                                     ints_dc[key4], rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = curve_interpolation.smooth_curve(tm_dc[key2], avg_pos_dc[key7],
                                                                                     ints_dc[key4], rep=1)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    line_styles = ['b--', 'r--', 'g--', 'm--']
    labels1 = []
    labels2 = []
    for i, value in enumerate(wsweep):
        label = 'width = %s' % value
        xname = 'smooth_tm' + str(i + 1)
        yname1 = 'smooth_pos' + str(i + 1)
        yname2 = 'smooth_avg' + str(i + 1)
        labels1.append(label)
        labels2.append(label)

        if value == 0.1:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label='width = %s (base)' % value)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label='width = %s (base)' % value)
        else:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label=label)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label=label)

    ax1.legend(loc='best', prop={'size': legend_size})
    # Change below to reflect title
    ax1.set_title("Varying width of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': legend_size})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')
    plt.tight_layout()
    if filename is not None:
        plt.savefig(f'{target_directory}{filename}_width_var.pdf', format='pdf')
    else:
        plt.savefig(f'{target_directory}width_var.{file_format}', format=file_format)
        print(f"Warning: filename_p is None. Saving as width_var.{file_format}")


def plot_traj_incvar(processed_adata, filename=None, isweep: List[float] = None, fig_width=5, fig_height=7,
                     legend_size=medium_size, file_format='pdf'):
    # Ensure that the target directory exists
    if filename is not None:
        target_directory = f'figures/{filename}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Warning: filename is None. Saving to 'figures/' directory.")

    # Default isweep value if not provided
    default_isweep = np.linspace(0.05, 0.2, num=4)
    default_isweep = np.around(default_isweep, 4)

    # Use provided isweep or default
    isweep = isweep if isweep is not None else default_isweep

    # Check if isweep is a vector with exactly 4 values
    if len(isweep) != 4:
        raise ValueError("isweep must be a vector with exactly 4 values.")

    # create dictionaries to store values for plotting
    max_eig_dc = {}
    tm_dc = {}
    pos_dc = {}
    ints_dc = {}
    npoints_dc = {}
    avg_max_dc = {}
    avg_pos_dc = {}
    std_max_dc = {}
    std_pos_dc = {}

    smooth_tm_dc = {}
    smooth_avg_dc = {}
    smooth_pos_dc = {}

    for i, value in enumerate(isweep):
        key1 = "max_eig" + str(i + 1)
        key2 = "tm" + str(i + 1)
        key3 = "pos" + str(i + 1)
        key4 = "ints" + str(i + 1)
        key5 = "npoints" + str(i + 1)
        key6 = "avg_max" + str(i + 1)
        key7 = "avg_pos" + str(i + 1)
        key8 = "std_max" + str(i + 1)
        key9 = "std_pos" + str(i + 1)
        key10 = "smooth_tm" + str(i + 1)
        key11 = "smooth_avg" + str(i + 1)
        key12 = "smooth_pos" + str(i + 1)

        max_eig_dc[key1] = processed_adata.uns['var_par']['inc=%s' % value]['largest']
        tm_dc[key2] = processed_adata.uns['var_par']['inc=%s' % value]['time']
        pos_dc[key3] = processed_adata.uns['var_par']['inc=%s' % value]['positive']
        ints_dc[key4] = processed_adata.uns['var_par']['inc=%s' % value]['pst_interval']
        npoints_dc[key5] = processed_adata.uns['var_par']['inc=%s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = curve_interpolation.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = curve_interpolation.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    line_styles = ['b--', 'r--', 'g--', 'm--']
    labels1 = []
    labels2 = []
    for i, value in enumerate(isweep):
        label = 'inc = %s' % value
        xname = 'smooth_tm' + str(i + 1)
        yname1 = 'smooth_pos' + str(i + 1)
        yname2 = 'smooth_avg' + str(i + 1)
        labels1.append(label)
        labels2.append(label)
        if value == 0.05:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label='inc = %s (base)' % value)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label='inc = %s (base)' % value)
        else:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label=label)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label=label)

    ax1.legend(loc='best', prop={'size': legend_size})
    # Change below to reflect title
    ax1.set_title("Varying increase of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': legend_size})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    if filename is not None:
        plt.savefig(f'{target_directory}{filename}_inc_var.{file_format}', format={file_format})
    else:
        plt.savefig(f'{target_directory}inc_var.{file_format}', format=file_format)
        print(f"Warning: filename_p is None. Saving as inc_var.{file_format}")


def plot_traj_nsimvar(adata, filename=None, nsweep: List[float] = None, fig_width=5, fig_height=7,
                      legend_size=medium_size, file_format='pdf'):
    # Ensure that the target directory exists
    if filename is not None:
        target_directory = f'figures/{filename}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Warning: filename is None. Saving to 'figures/' directory.")

    # Default nsweep value if not provided
    default_nsweep = np.linspace(10, 40, num=4)
    default_nsweep = np.around(default_nsweep, 0)

    # Use provided nsweep or default
    nsweep = nsweep if nsweep is not None else default_nsweep

    # Check if nsweep is a vector with exactly 4 values
    if len(nsweep) != 4:
        raise ValueError("nsweep must be a vector with exactly 4 values.")

    # create dictionaries to store values for plotting
    max_eig_dc = {}
    tm_dc = {}
    pos_dc = {}
    ints_dc = {}
    npoints_dc = {}
    avg_max_dc = {}
    avg_pos_dc = {}
    std_max_dc = {}
    std_pos_dc = {}

    smooth_tm_dc = {}
    smooth_avg_dc = {}
    smooth_pos_dc = {}

    for i, value in enumerate(nsweep):
        key1 = "max_eig" + str(i + 1)
        key2 = "tm" + str(i + 1)
        key3 = "pos" + str(i + 1)
        key4 = "ints" + str(i + 1)
        key5 = "npoints" + str(i + 1)
        key6 = "avg_max" + str(i + 1)
        key7 = "avg_pos" + str(i + 1)
        key8 = "std_max" + str(i + 1)
        key9 = "std_pos" + str(i + 1)
        key10 = "smooth_tm" + str(i + 1)
        key11 = "smooth_avg" + str(i + 1)
        key12 = "smooth_pos" + str(i + 1)

        max_eig_dc[key1] = adata.uns['var_par']['nsim=%s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['nsim=%s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['nsim=%s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['nsim=%s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['nsim=%s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = curve_interpolation.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = curve_interpolation.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    line_styles = ['b--', 'r--', 'g--', 'm--']
    labels1 = []
    labels2 = []
    for i, value in enumerate(nsweep):
        label = 'nsim = %s' % value
        xname = 'smooth_tm' + str(i + 1)
        yname1 = 'smooth_pos' + str(i + 1)
        yname2 = 'smooth_avg' + str(i + 1)
        labels1.append(label)
        labels2.append(label)

        if value == 10:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label='nsim = %s (base)' % value)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label='nsim = %s (base)' % value)
        else:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label=label)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label=label)

    ax1.legend(loc='best', prop={'size': legend_size})
    # Change below to reflect title
    ax1.set_title("Varying nsim of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': legend_size})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    if filename is not None:
        plt.savefig(f'{target_directory}{filename}_nsim_var.{file_format}', format={file_format})
    else:
        plt.savefig(f'{target_directory}nsim_var.{file_format}', format=file_format)
        print(f"Warning: filename_p is None. Saving as nsim_var.{file_format}")


def plot_traj_fracvar(adata, filename=None, fsweep: List[float] = None, fig_width=5, fig_height=7,
                      legend_size=medium_size, file_format='pdf'):
    # Ensure that the target directory exists
    if filename is not None:
        target_directory = f'figures/{filename}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Warning: filename is None. Saving to 'figures/' directory.")

    # Default fsweep value if not provided
    default_fsweep = np.linspace(0.3, 0.9, num=4)
    default_fsweep = np.around(default_fsweep, 4)

    # Use provided fsweep or default
    fsweep = fsweep if fsweep is not None else default_fsweep

    # Check if fsweep is a vector with exactly 4 values
    if len(fsweep) != 4:
        raise ValueError("fsweep must be a vector with exactly 4 values.")

    # create dictionaries to store values for plotting
    max_eig_dc = {}
    tm_dc = {}
    pos_dc = {}
    ints_dc = {}
    npoints_dc = {}
    avg_max_dc = {}
    avg_pos_dc = {}
    std_max_dc = {}
    std_pos_dc = {}

    smooth_tm_dc = {}
    smooth_avg_dc = {}
    smooth_pos_dc = {}

    for i, value in enumerate(fsweep):
        key1 = "max_eig" + str(i + 1)
        key2 = "tm" + str(i + 1)
        key3 = "pos" + str(i + 1)
        key4 = "ints" + str(i + 1)
        key5 = "npoints" + str(i + 1)
        key6 = "avg_max" + str(i + 1)
        key7 = "avg_pos" + str(i + 1)
        key8 = "std_max" + str(i + 1)
        key9 = "std_pos" + str(i + 1)
        key10 = "smooth_tm" + str(i + 1)
        key11 = "smooth_avg" + str(i + 1)
        key12 = "smooth_pos" + str(i + 1)

        max_eig_dc[key1] = adata.uns['var_par']['frac=%s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['frac=%s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['frac=%s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['frac=%s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['frac=%s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = curve_interpolation.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = curve_interpolation.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    line_styles = ['b--', 'r--', 'g--', 'm--']
    labels1 = []
    labels2 = []
    for i, value in enumerate(fsweep):
        label = 'frac = %s' % value
        xname = 'smooth_tm' + str(i + 1)
        yname1 = 'smooth_pos' + str(i + 1)
        yname2 = 'smooth_avg' + str(i + 1)
        labels1.append(label)
        labels2.append(label)

        if value == 0.9:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label='frac = %s (base)' % value)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label='frac = %s (base)' % value)
        else:
            ax1.plot(smooth_tm_dc[xname], smooth_pos_dc[yname1], line_styles[i], label=label)
            ax2.plot(smooth_tm_dc[xname], smooth_avg_dc[yname2], line_styles[i], label=label)

    ax1.legend(loc='best', prop={'size': legend_size})
    # Change below to reflect title
    ax1.set_title("Varying frac of cells in window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': legend_size})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()

    if filename is not None:
        plt.savefig(f'{target_directory}{filename}_frac_var.{file_format}', format={file_format})
    else:
        plt.savefig(f'{target_directory}frac_var.{file_format}', format=file_format)
        print(f"Warning: filename_p is None. Saving as frac_var.{file_format}")
