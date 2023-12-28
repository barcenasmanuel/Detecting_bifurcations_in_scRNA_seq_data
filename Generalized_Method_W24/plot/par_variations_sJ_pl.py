import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import scanpy as sc

small_size = 6
medium_size = 8
large_size = 12


def plot_traj_widthvar(adata, wsweep):
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

        max_eig_dc[key1] = adata.uns['var_par']['width= %s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['width= %s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['width= %s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['width= %s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['width= %s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = tools.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = tools.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(5, 7))
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

    ax1.legend(loc='best', prop={'size': 8})
    # Change below to reflect title
    ax1.set_title("Varying width of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': 8})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    plt.savefig('figures/width_var.pdf', format='pdf')


def plot_traj_incvar(adata, isweep):
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

        max_eig_dc[key1] = adata.uns['var_par']['inc= %s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['inc= %s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['inc= %s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['inc= %s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['inc= %s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = tools.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = tools.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(5, 7))
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

    ax1.legend(loc='best', prop={'size': 8})
    # Change below to reflect title
    ax1.set_title("Varying increase of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': 8})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    plt.savefig('figures/inc_var.pdf', format='pdf')


def plot_traj_fracvar(adata, fsweep):
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

        max_eig_dc[key1] = adata.uns['var_par']['frac= %s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['frac= %s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['frac= %s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['frac= %s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['frac= %s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = tools.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = tools.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(5, 7))
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

    ax1.legend(loc='best', prop={'size': 8})
    # Change below to reflect title
    ax1.set_title("Varying frac of cells in window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': 8})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    plt.savefig('figures/frac_var.pdf', format='pdf')


def plot_traj_nsimvar(adata, nssweep):
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

    for i, value in enumerate(nssweep):
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

        max_eig_dc[key1] = adata.uns['var_par']['nsim= %s' % value]['largest']
        tm_dc[key2] = adata.uns['var_par']['nsim= %s' % value]['time']
        pos_dc[key3] = adata.uns['var_par']['nsim= %s' % value]['positive']
        ints_dc[key4] = adata.uns['var_par']['nsim= %s' % value]['pst_interval']
        npoints_dc[key5] = adata.uns['var_par']['nsim= %s' % value]['npoints']

        avg_max_dc[key6], avg_pos_dc[key7] = np.mean(max_eig_dc[key1], axis=1), np.mean(pos_dc[key3], axis=1)
        std_max_dc[key8], std_pos_dc[key9] = np.std(max_eig_dc[key1], axis=1), np.std(pos_dc[key3], axis=1)

        smooth_tm_dc[key10], smooth_avg_dc[key11] = tools.smooth_curve(tm_dc[key2], avg_max_dc[key6], ints_dc[key4],
                                                                       rep=1)
        smooth_tm_dc[key10], smooth_pos_dc[key12] = tools.smooth_curve(tm_dc[key2], avg_pos_dc[key7], ints_dc[key4],
                                                                       rep=1)

    fig = plt.figure(figsize=(5, 7))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    line_styles = ['b--', 'r--', 'g--', 'm--']
    labels1 = []
    labels2 = []
    for i, value in enumerate(nssweep):
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

    ax1.legend(loc='best', prop={'size': 8})
    # Change below to reflect title
    ax1.set_title("Varying nsim of window on spliceJAC")

    ax1.set_ylabel('Positive eigenvalues')
    ax1.set_xticks([])

    ax2.legend(loc='best', prop={'size': 8})
    ax2.set_ylabel('Largest eigenvalue')
    ax2.set_xlabel('Pseudotime')

    plt.tight_layout()
    plt.savefig('figures/nsim_var.pdf', format='pdf')
