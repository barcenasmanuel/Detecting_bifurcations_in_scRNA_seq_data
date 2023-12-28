import matplotlib
import matplotlib.pyplot as plt
import numpy as np

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


def plot_comm(adata, num_com_list):
    tm = adata.uns['Jacobian']['time']

    matplotlib.rc('xtick', labelsize=small_size)
    matplotlib.rc('ytick', labelsize=small_size)
    fig = plt.figure(figsize=(3.54, 3.54), dpi=600)  # original

    plt.plot(tm, num_com_list, 'mo-')  # bo- #b--
    plt.ylabel('Number of Communities', fontsize=medium_size)
    plt.xlabel('Pseudotime', fontsize=medium_size)

    num_ticks = 6
    xtickpos = np.linspace(min(tm), max(tm), num_ticks)
    ytickpos = np.linspace(0, max(num_com_list) + 3, num_ticks + 1)
    plt.xticks(xtickpos)
    plt.yticks(ytickpos)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%d'))  # set y as integers
    # OVCA420_TGFB1
    # No Gene expression
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_GirvanNewman2.pdf', format='pdf')
    # Gene expression
    # girvan newman
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_GirvanNewman1.pdf', format='pdf', bbox_inches='tight')
    # greedy_modularity
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_greedyModularity1.pdf', format='pdf', bbox_inches='tight')
    # OVCA420_EGF
    # plt.savefig('figures/1st_draft/OVCA420_EGF/Ncomm_vs_pst_GirvanNewman1_OVCA420_EGF.pdf', format='pdf', bbox_inches='tight')
    # plt.savefig('figures/1st_draft/OVCA420_EGF/Ncomm_vs_pst_greedyModularity1_OVCA420_EGF.pdf', format='pdf', bbox_inches='tight')

def plot_comm_multi(adata, num_com_list1, num_com_list2, filename_p):
    tm = adata.uns['Jacobian']['time']

    matplotlib.rc('xtick', labelsize=small_size)
    matplotlib.rc('ytick', labelsize=small_size)
    fig = plt.figure(figsize=(3.54, 3.54))  # original #dpi=600

    plt.plot(tm, num_com_list1, 'co-', label='Girvan_Newman')  #
    plt.plot(tm, num_com_list2, 'mo-', label='Greedy_modularity')  # bo- #b--
    plt.ylabel('Number of Communities', fontsize=8)
    plt.xlabel('Pseudotime', fontsize=8)
    plt.legend(loc='best', prop={'size': 6})

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
    # OVCA420_TGFB1
    # No Gene expression
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_GirvanNewman2.pdf', format='pdf')
    # Gene expression
    # girvan newman
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_GirvanNewman1.pdf', format='pdf', bbox_inches='tight')
    # greedy_modularity
    # plt.savefig('figures/1st_draft/Ncomm_vs_pst_greedyModularity1.pdf', format='pdf', bbox_inches='tight')
    # Multiple
    plt.savefig(f'figures/1st_draft/{filename_p}/Ncomm_vs_pst_multi_methods2_{filename_p}.png', format='png',
                bbox_inches='tight', dpi=300)
    # greedy_modularity
    # OVCA420_EGF
    # plt.savefig('figures/1st_draft/OVCA420_EGF/Ncomm_vs_pst_GirvanNewman1_OVCA420_EGF.pdf', format='pdf', bbox_inches='tight')
    # plt.savefig('figures/1st_draft/OVCA420_EGF/Ncomm_vs_pst_greedyModularity1_OVCA420_EGF.pdf', format='pdf', bbox_inches='tight')
    # Multiple
    # plt.savefig('figures/1st_draft/OVCA420_EGF/Ncomm_vs_pst_multi_methods1_OVCA420_EGF.pdf', format='pdf', bbox_inches='tight')


def plot_distribution(G_list, filename_p):
    all_edge_weights = []
    print(f'Running comm_plots.plot_distribution')
    # Iterate over the GRNs and plot their edge weight distributions with KDE curves
    for i, G in enumerate(G_list):
        if i % 3 == 0:
            edge_weights = [data['weight'] for _, _, data in G.edges(data=True) if 'weight' in data]
            all_edge_weights.extend(edge_weights)

    min_edge_weight = min(all_edge_weights)
    max_edge_weight = max(all_edge_weights)

    matplotlib.rc('xtick', labelsize=small_size)
    matplotlib.rc('ytick', labelsize=small_size)
    plt.figure(figsize=(3.54, 3.54))

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

    plt.xlabel("Edge Weight", fontsize=medium_size)
    plt.ylabel("Density", fontsize=medium_size)
    plt.grid(True)
    plt.legend(fontsize=4)
    # OVCA420_TGFB1
    plt.xlim([-0.75, 0.75])  # for geneweight expression OVCA420_TGFB1
    # plt.xlim([-0.25,0.25]) #for weight expression OVCA420_TGFB1
    # OVCA420_EGF
    # plt.xlim([-0.21,0.21])
    # plt.xlim([-0.1,0.1])
    # OVCA420_TGFB1
    # plt.savefig('figures/1st_draft/geneWeight_distribution.pdf', format='pdf', bbox_inches='tight')
    # plt.savefig('figures/1st_draft/Weight_distribution.pdf', format='pdf', bbox_inches='tight')
    # OVCA420_EGF
    # plt.savefig('figures/1st_draft/OVCA420_EGF/geneWeight_distribution_OVCA420_EGF_A.pdf', format='pdf', bbox_inches='tight')
    plt.savefig(f'figures/1st_draft/{filename_p}/Weight_distribution_{filename_p}.png', format='png',
                bbox_inches='tight', dpi=300)

    print(f'End comm_plots.plot_distribution \nsee figure folder for figure')
