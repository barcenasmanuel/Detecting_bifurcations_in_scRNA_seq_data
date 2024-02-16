import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import networkx as nx
from scipy.stats import gaussian_kde

small_size = 6
medium_size = 8
large_size = 12


def comm_nodes(G):
    communities = nx.community.asyn_lpa_communities(G, seed=1)
    node_groups = []
    '''
    Create the seed to fix position of Communities generated 
    '''
    seed_value = 7
    pos_nr = nx.spring_layout(G, seed=seed_value)
    '''
    for com in communities:
        node_groups.append(list(com)) # Convert the community to a list and append it to node_groups
    '''
    node_groups = [list(community) for community in communities]
    n = len(node_groups)
    # color_map = cm.rainbow(np.linspace(0,1,n))
    color_map = plt.cm.tab20.colors[:n]

    node_colors = {}
    for i, community in enumerate(node_groups):
        for node in community:
            node_colors[node] = color_map[i]

    node_colors_list = [node_colors[node] for node in G.nodes]
    # draw_networkx appears to randomize the pos (position)
    # nx.draw_networkx(G, node_color=node_colors_list, with_labels=True, font_size= 6) #font_size=8
    nx.draw_networkx(G, pos=pos_nr, node_color=node_colors_list, with_labels=True, font_size=4)  # font_size=8

    plt.savefig('figures/nodes_draw_test4c.pdf', format='pdf')


def plot_comm(adata, num_com_list, filename_p=None, plot_tick=small_size, axis_label_size=medium_size,
              fig_width=3.54, fig_height=3.54, fig_dpi=100, method_key='girvan_newman',
              grn_type='gene_exp_weights', file_format='png'):
    """
    Plot the number of communities over pseudotime.

    Parameters:
    - adata (AnnData): Annotated data object containing pseudotime information.
    - num_com_list (list): List of the number of communities corresponding to each pseudotime point.
    - filename_p (str, optional): Output file name and directory (without extension).
    - plot_tick (int): Font size for tick labels.
    - axis_label_size (int): Font size for axis labels.
    - fig_width (float): Width of the figure in inches.
    - fig_height (float): Height of the figure in inches.
    - fig_dpi (int): Dots per inch for the figure resolution.
    - method_key (str): Key representing the community detection method.
    - grn_type (str): Type of gene regulatory network (GRN) weights used.
    - file_format (str): Output file format for saving the plot. (i.e: png, pdf, etc)

    Returns:
    None
    """
    print(f'Running plot_comm')
    # Ensure that the target directory exists
    if filename_p is not None:
        target_directory = f'figures/{filename_p}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Note: filename_p is None. Saving to 'figures/' directory.")

    tm = adata.uns['Jacobian']['time']

    matplotlib.rc('xtick', labelsize=plot_tick)
    matplotlib.rc('ytick', labelsize=plot_tick)
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=fig_dpi)  # original

    plt.plot(tm, num_com_list, 'mo-')  # bo- #b--
    plt.ylabel('Number of Communities', fontsize=axis_label_size)
    plt.xlabel('Pseudotime', fontsize=axis_label_size)

    num_ticks = 6
    xtickpos = np.linspace(min(tm), max(tm), num_ticks)
    ytickpos = np.linspace(0, max(num_com_list) + 3, num_ticks + 1)
    plt.xticks(xtickpos)
    plt.yticks(ytickpos)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%d'))  # set y as integers

    if grn_type is not None:
        plt.savefig(f'{target_directory}{method_key}_{filename_p}_{grn_type}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)
    else:
        plt.savefig(f'{target_directory}{method_key}_{filename_p}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)


def plot_comm_multi(adata, num_com_list1, num_com_list2, filename_p, plot_tick=small_size, axis_label_size=medium_size,
                    legend_size=small_size, fig_width=3.54, fig_height=3.54, fig_dpi=100,
                    method_key1='girvan_newman', method_key2='greedy_modularity', grn_type=None, file_format='pdf'):
    """
        Plot the number of communities over pseudotime for two different community detection methods.

        Parameters:
        - adata (AnnData): Annotated data object containing pseudotime information.
        - num_com_list1 (list): List of the number of communities corresponding to method_key1.
        - num_com_list2 (list): List of the number of communities corresponding to method_key2.
        - filename_p (str): Output file name (without extension) and directory.
        - plot_tick (int): Font size for tick labels.
        - axis_label_size (int): Font size for axis labels.
        - legend_size (int): Font size for legend.
        - fig_width (float): Width of the figure in inches.
        - fig_height (float): Height of the figure in inches.
        - fig_dpi (int): Dots per inch for the figure resolution.
        - method_key1 (str): Key representing the first community detection method.
        - method_key2 (str): Key representing the second community detection method.
        - grn_type (str, optional): Type of gene regulatory network (GRN) weights used. Default is None.
        - file_format (str): Output file format for saving the plot. Default is 'pdf'.

        Returns:
        None

        Saves a plot showing the number of communities over pseudotime for two different methods.
        The plot is saved in the specified directory with a filename based on method keys and GRN type.

        If filename_p is None, the plot is saved to the 'figures/' directory.
        If grn_type is provided, it is included in the filename for better identification.
    """
    # Ensure that the target directory exists
    if filename_p is not None:
        target_directory = f'figures/{filename_p}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Note: filename_p is None. Saving to 'figures/' directory.")

    tm = adata.uns['Jacobian']['time']

    matplotlib.rc('xtick', labelsize=plot_tick)
    matplotlib.rc('ytick', labelsize=plot_tick)
    fig = plt.figure(figsize=(fig_width, fig_height))  # original #dpi=600

    plt.plot(tm, num_com_list1, 'co-', label=method_key1)  #
    plt.plot(tm, num_com_list2, 'mo-', label=method_key2)  # bo- #b--
    plt.ylabel('Number of Communities', fontsize=axis_label_size)
    plt.xlabel('Pseudotime', fontsize=axis_label_size)
    plt.legend(loc='best', prop={'size': legend_size})

    print(f'num_com_list1: {num_com_list1}')
    print(f'num_com_list2: {num_com_list2}')

    max1 = max(num_com_list1)
    max2 = max(num_com_list2)
    if max1 < max2:
        xmax = max1
    else:
        xmax = max2

    num_ticks = 6
    xtickpos = np.linspace(min(tm), max(tm), num_ticks)
    ytickpos = np.linspace(0, xmax + 3, num_ticks + 1)
    plt.xticks(xtickpos)
    plt.yticks(ytickpos)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%d'))  # set y as integers

    # Multiple
    if grn_type is not None:
        plt.savefig(f'{target_directory}Ncomm_vs_pst_multi_methods_{filename_p}_{grn_type}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)
    else:
        plt.savefig(f'{target_directory}Ncomm_vs_pst_multi_methods_{filename_p}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)


def plot_distribution(G_list, filename_p, plot_tick=small_size, axis_label_size=medium_size,legend_size=4,
                      fig_width=3.54, fig_height=3.54, fig_dpi=100, grn_type=None, xlim_p=None, file_format='png'):
    """
    Plot the distribution of edge weights for a list of gene regulatory networks (GRNs).

    Parameters:
        G_list (list): A list of networkx graphs representing GRNs.
        filename_p (str): Output file name and directory (without extension).
        plot_tick (int): Font size for tick labels.
        axis_label_size (int): Font size for axis labels.
        legend_size (int): Font size for legend.
        fig_width (float): Width of the figure in inches.
        fig_height (float): Height of the figure in inches.
        fig_dpi (int): Dots per inch for the figure resolution.
        grn_type (str): Type of GRN, if applicable.
        xlim_p (tuple): Tuple specifying the range for the x-axis limit.
        file_format (str): Output file format for saving the plot.

    Returns:
        None
    """
    # Ensure that the target directory exists
    if filename_p is not None:
        target_directory = f'figures/{filename_p}/'
        os.makedirs(target_directory, exist_ok=True)
    else:
        target_directory = 'figures/'
        print(f"Note: filename_p is None. Saving to 'figures/' directory.")

    all_edge_weights = []
    print(f'Running plot_distribution')
    # Iterate over the GRNs and plot their edge weight distributions with KDE curves
    for i, G in enumerate(G_list):
        if i % 3 == 0:
            edge_weights = [data['weight'] for _, _, data in G.edges(data=True) if 'weight' in data]
            all_edge_weights.extend(edge_weights)

    min_edge_weight = min(all_edge_weights)
    max_edge_weight = max(all_edge_weights)

    matplotlib.rc('xtick', labelsize=plot_tick)
    matplotlib.rc('ytick', labelsize=plot_tick)
    plt.figure(figsize=(fig_width, fig_height), dpi=fig_dpi)

    for i, G in enumerate(G_list):
        if i % 3 == 0:
            edge_weights = [data['weight'] for _, _, data in G.edges(data=True) if 'weight' in data]
            # Calculate KDE curve
            kde = gaussian_kde(edge_weights)  # , bw_method = common_bandwidth [add this to restructure distributions]
            x_vals = np.linspace(min_edge_weight, max_edge_weight, 100)
            y_vals = kde(x_vals)

            # Normalize the KDE curve
            y_vals /= np.trapz(y_vals, x_vals)

            # Plot the KDE curve
            plt.plot(x_vals, y_vals, lw=2, label=f'KDE GRN {i}')

    plt.xlabel("Edge Weight", fontsize=axis_label_size)
    plt.ylabel("Density", fontsize=axis_label_size)
    plt.grid(True)
    plt.legend(fontsize=legend_size)
    if xlim_p is not None:
        plt.xlim(xlim_p)

    # OVCA420_TGFB1
    # plt.xlim([-0.75, 0.75])  # for geneweight expression OVCA420_TGFB1
    # plt.xlim([-0.25,0.25]) #for weight expression OVCA420_TGFB1
    # OVCA420_EGF
    # plt.xlim([-0.21,0.21])
    # plt.xlim([-0.1,0.1])
    # OVCA420_TGFB1
    # plt.savefig('figures/1st_draft/geneWeight_distribution.pdf', format='pdf', bbox_inches='tight')
    # plt.savefig('figures/1st_draft/Weight_distribution.pdf', format='pdf', bbox_inches='tight')
    # OVCA420_EGF
    # plt.savefig('figures/1st_draft/OVCA420_EGF/geneWeight_distribution_OVCA420_EGF_A.pdf', format='pdf', bbox_inches='tight')

    if grn_type is not None:
        plt.savefig(f'{target_directory}Weight_distribution_{filename_p}_{grn_type}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)
    else:
        plt.savefig(f'{target_directory}Weight_distribution_{filename_p}.{file_format}',
                    format=file_format, bbox_inches='tight', dpi=fig_dpi)

    plt.savefig(f'{target_directory}Weight_distribution_{filename_p}.{file_format}', format=file_format,
                bbox_inches='tight', dpi=fig_dpi)
