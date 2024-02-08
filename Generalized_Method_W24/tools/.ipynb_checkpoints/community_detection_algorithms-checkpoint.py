import numpy as np
import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
from networkx import edge_betweenness_centrality as betweenness
from operator import itemgetter


def calc_grn(adata, genes=None, index=0, weight_quantile=0.9):
    """
        Calculate a Gene Regulatory Network (GRN) from the Jacobian matrix calculated from grn_inference function.

        Parameters:
            adata (AnnData): Annotated data matrix.
            genes (list, optional): List of gene names. If not provided, it uses all genes in `adata`.
            index (int, optional): Index of the Jacobian matrix.
            weight_quantile (float, optional): Weight quantile for thresholding edges.

        Returns:
            networkx.DiGraph: Directed graph representing the GRN.

    """
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
    """
        Calculate a Gene Regulatory Network (GRN) with mean gene expression multiplied to edge weights.

        Parameters:
            adata (AnnData): Annotated data matrix.
            genes (list, optional): List of gene names. If not provided, it uses all genes in `adata`.
            index (int, optional): Index of the Jacobian matrix.
            weight_quantile (float, optional): Weight quantile for thresholding edges.

        Returns:
            networkx.DiGraph: Directed graph representing the GRN with gene expression weights.

    """
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


def G_listgen_a(adata, jac_list, wq=0.9):
    """
        Generate a list of Gene Regulatory Networks (GRNs) using the calc_grn function.

        Parameters:
            adata (AnnData): Annotated data matrix.
            jac_list (list): List of Jacobian matrices.
            wq (float, optional): Weight quantile for thresholding edges.

        Returns:
            list: List of GRNs.

    """
    G_list = []

    for i, jac_value in enumerate(jac_list):
        G = calc_grn(adata, index=i, weight_quantile=wq)  # wq=0.9
        G_list.append(G)

    return G_list


def G_listgen_b(adata, jac_list, wq=0.9):
    """
    Generate a list of Gene Regulatory Networks (GRNs) with gene expression weights using the grn_geneweight_exp function.

    Parameters:
        adata (AnnData): Annotated data matrix.
        jac_list (list): List of Jacobian matrices.
        wq (float, optional): Weight quantile for thresholding edges.

    Returns:
        list: List of GRNs with gene expression weights.

    """
    G_list = []

    for i, jac_value in enumerate(jac_list):
        G = grn_geneweight_exp(adata, index=i, weight_quantile=wq)  # wq=0.9
        G_list.append(G)

    return G_list


def between_method(G):
    """
        Calculate the node with the highest betweenness centrality in a graph.

        Parameters:
            G (networkx.Graph): Graph for calculating betweenness centrality.

        Returns:
            hashable: Node with the highest betweenness centrality.

    """
    centrality = betweenness(G, weight='weight')
    return max(centrality, key=centrality.get)


def cd_g_w(adata, jac_list=None, g_list_method="G_list_gen_b", grn_list=None,
           most_valuable_edge=between_method, wq=0.9):
    """
        Detect communities in a graph based on weights of graph using girvan_newman method.

        Parameters:
        - adata: AnnData object
        - jac_list: List of Jacobian matrices
        - g_list_method: Method to generate the graph list, choose from ["G_list_gen_b", "G_list_gen_a", "Other"]
        - grn_list: Graph list for "Other" method
        - most_valuable_edge: Edge metric for community detection using betweenness centrality method
        - wq: Weight quantile for graph generation

        Returns:
        - Community_list: List of communities detected for each graph
        - Num_com_list: List of the number of communities for each graph
    """
    if jac_list is None:
        jac_list = adata.uns['Jacobian']['jacobians']

    # Parameter validation
    if g_list_method not in ["G_list_gen_b", "G_list_gen_a", "Other"]:
        raise ValueError("Invalid value for g_list_method")

    if grn_list is None and g_list_method == "Other":
        raise ValueError("Please provide a grn_list when using 'Other' as the g_list_method.")

    # Graph list generation
    if g_list_method == "G_list_gen_b":
        G_list = G_listgen_b(adata, jac_list, wq)
    elif g_list_method == "G_list_gen_a":
        G_list = G_listgen_a(adata, jac_list, wq)
    elif g_list_method == "Other":
        G_list = grn_list

    Community_list = []
    Num_com_list = []

    # Community detection loop
    for i, g_list in enumerate(G_list):
        comm_gen = girvan_newman(g_list, most_valuable_edge=most_valuable_edge)
        communities = [list(community) for community in next(comm_gen)]
        Community_list.append(communities)

        num_communities = len(communities)
        Num_com_list.append(num_communities)

        community_sizes = [len(community) for community in communities]

        average_size = sum(community_sizes) / len(community_sizes)

    return Community_list, Num_com_list


def cd_g_nw(adata, jac_list=None, g_list_method="G_list_gen_b", grn_list=None, wq=0.9):
    """
            Detect communities in a graph using girvan_newman method with all edges having the same weight.

            Parameters:
            - adata: AnnData object
            - jac_list: List of Jacobian matrices
            - g_list_method: Method to generate the graph list, choose from ["G_list_gen_b", "G_list_gen_a", "Other"]
            - grn_list: Graph list for "Other" method
            - wq: Weight quantile for graph generation

            Returns:
            - Community_list: List of communities detected for each graph
            - Num_com_list: List of the number of communities for each graph
        """
    if jac_list is None:
        jac_list = adata.uns['Jacobian']['jacobians']

    # Parameter validation
    if g_list_method not in ["G_list_gen_b", "G_list_gen_a", "Other"]:
        raise ValueError("Invalid value for g_list_method")

    if grn_list is None and g_list_method == "Other":
        raise ValueError("Please provide a grn_list when using 'Other' as the g_list_method.")

    # Graph list generation
    if g_list_method == "G_list_gen_b":
        G_list = G_listgen_b(adata, jac_list, wq)
    elif g_list_method == "G_list_gen_a":
        G_list = G_listgen_a(adata, jac_list, wq)
    elif g_list_method == "Other":
        G_list = grn_list

    Community_list = []
    Num_com_list = []

    # Community detection loop
    for i, g_list in enumerate(G_list):
        comm_gen = girvan_newman(g_list)
        communities = [list(community) for community in next(comm_gen)]
        Community_list.append(communities)

        num_communities = len(communities)
        Num_com_list.append(num_communities)

        community_sizes = [len(community) for community in communities]

        average_size = sum(community_sizes) / len(community_sizes)

    return Community_list, Num_com_list


def cd_grM_w(adata, jac_list=None, g_list_method="G_list_gen_b", grn_list=None, wq=0.9):
    """
        Detect communities in gene regulatory networks (GRNs) using weighted greedy modularity optimization.

        Parameters:
            adata (AnnData): Annotated data matrix.
            jac_list (list, optional): List of Jacobian matrices. If not provided, it uses the Jacobians stored in `adata`.
            g_list_method (str, optional): Method for generating the graph list. Choices: "G_list_gen_b", "G_list_gen_a", "Other".
            grn_list (list, optional): List of graphs when using "Other" as the graph list generation method.
            wq (float, optional): Weight quantile for generating weighted graphs.

        Returns:
            tuple: Tuple containing two lists: `Community_list` (communities detected for each GRN) and `Num_com_list`
            (number of communities detected for each GRN).

        Raises:
            ValueError: If an invalid value is provided for `g_list_method` or if `grn_list` is not provided when using "Other".

    """
    if jac_list is None:
        jac_list = adata.uns['Jacobian']['jacobians']

    # Parameter validation
    if g_list_method not in ["G_list_gen_b", "G_list_gen_a", "Other"]:
        raise ValueError("Invalid value for g_list_method")

    if grn_list is None and g_list_method == "Other":
        raise ValueError("Please provide a grn_list when using 'Other' as the g_list_method.")

    # Graph list generation
    if g_list_method == "G_list_gen_b":
        G_list = G_listgen_b(adata, jac_list, wq)
    elif g_list_method == "G_list_gen_a":
        G_list = G_listgen_a(adata, jac_list, wq)
    elif g_list_method == "Other":
        G_list = grn_list

    Community_list = []
    Num_com_list = []

    # Community detection loop
    for i, g_list in enumerate(G_list):
        comm_gen = nx.community.greedy_modularity_communities(g_list, weight='weight')  # with weights
        communities = [list(community) for community in comm_gen]
        Community_list.append(communities)

        num_communities = len(communities)
        Num_com_list.append(num_communities)

        community_sizes = [len(community) for community in communities]

        average_size = sum(community_sizes) / len(community_sizes)

    return Community_list, Num_com_list


def cd_grM_nw(adata, jac_list=None, g_list_method="G_list_gen_b", grn_list=None, wq=0.9):
    """
        Detect communities in gene regulatory networks (GRNs) using unweighted modularity optimization.

        Parameters:
            adata (AnnData): Annotated data matrix.
            jac_list (list, optional): List of Jacobian matrices. If not provided, it uses the Jacobians stored in `adata`.
            g_list_method (str, optional): Method for generating the graph list. Choices: "G_list_gen_b", "G_list_gen_a", "Other".
            grn_list (list, optional): List of graphs when using "Other" as the graph list generation method.
            wq (float, optional): Weight quantile for generating weighted graphs.

        Returns:
            tuple: Tuple containing two lists: `Community_list` (communities detected for each GRN) and `Num_com_list`
            (number of communities detected for each GRN).

        Raises:
            ValueError: If an invalid value is provided for `g_list_method` or if `grn_list` is not provided when using "Other".

     """
    if jac_list is None:
        jac_list = adata.uns['Jacobian']['jacobians']

    # Parameter validation
    if g_list_method not in ["G_list_gen_b", "G_list_gen_a", "Other"]:
        raise ValueError("Invalid value for g_list_method")

    if grn_list is None and g_list_method == "Other":
        raise ValueError("Please provide a grn_list when using 'Other' as the g_list_method.")

    # Graph list generation
    if g_list_method == "G_list_gen_b":
        G_list = G_listgen_b(adata, jac_list, wq)
    elif g_list_method == "G_list_gen_a":
        G_list = G_listgen_a(adata, jac_list, wq)
    elif g_list_method == "Other":
        G_list = grn_list

    Community_list = []
    Num_com_list = []

    # Community detection loop
    for i, g_list in enumerate(G_list):
        comm_gen = nx.community.greedy_modularity_communities(g_list)  # with weights
        communities = [list(community) for community in comm_gen]
        Community_list.append(communities)

        num_communities = len(communities)
        Num_com_list.append(num_communities)

        community_sizes = [len(community) for community in communities]

        average_size = sum(community_sizes) / len(community_sizes)

    return Community_list, Num_com_list
