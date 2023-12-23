import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
from networkx import edge_betweenness_centrality as betweenness

def between_method(G):
    '''
    betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    #can only use a seed when k=value
    :param G: Is the GRN in window from spliceJAC
    :return:
    '''
    centrality = betweenness(G, weight='weight')
    return max(centrality, key=centrality.get)


def cd_g_w(G_list, most_valuable_edge=between_method):
    '''

    :param G_list:
    :param most_valuable_edge:
    :return:
    '''
    Community_list = []
    Num_com_list = []
    print(f'Running community_detection_algorithms.community_detection_girvannewman_weights')
    for i, g_list in enumerate(G_list):
        # girvan_newman(G, most_valuable_edge=None)
        comm_gen = girvan_newman(g_list, most_valuable_edge=most_valuable_edge)  # with weights; runs slower

        communities = [list(community) for community in next(comm_gen)]
        Community_list.append(communities)
        # print(f'Done with iteration:{i}',' \n','Communities detected:', '\n', communities)

    # print(f'-'*120, f'\n Communities:', Community_list)
    for i, c_val in enumerate(Community_list):
        print(f'Communities for jacobian_{i}: ', '\n', c_val, '\n')

        # Compute the number of communities
        num_communities = len(c_val)
        Num_com_list.append(num_communities)
        print(f'Number of communities: {num_communities}')

        # Compute the size of each community in c_val
        community_sizes = [len(community) for community in c_val]
        print(f'Sizes of communities: {community_sizes}')

        # Compute the average size of communities in c_val
        average_size = sum(community_sizes) / len(community_sizes)
        print(f'Average community size: {average_size}\n')
    print(f'End comm_tools.community_detection_girvannewman_weights')
    return Community_list, Num_com_list


def cd_g_nw(G_list):
    Community_list = []
    Num_com_list = []
    print(f'Running comm_tools.community_detection_girvannewman_noweights')
    for i, g_list in enumerate(G_list):
        # girvan_newman(G, most_valuable_edge=None)
        comm_gen = girvan_newman(G=g_list)  # no weights; runs faster

        communities = [list(community) for community in next(comm_gen)]
        Community_list.append(communities)
        # print(f'Done with iteration:{i}',' \n','Communities detected:', '\n', communities)

    # print(f'-'*120, f'\n Communities:', Community_list)
    for i, c_val in enumerate(Community_list):
        print(f'Communities for jacobian_{i}: ', '\n', c_val, '\n')

        # Compute the number of communities
        num_communities = len(c_val)
        Num_com_list.append(num_communities)
        print(f'Number of communities: {num_communities}')

        # Compute the size of each community in c_val
        community_sizes = [len(community) for community in c_val]
        print(f'Sizes of communities: {community_sizes}')

        # Compute the average size of communities in c_val
        average_size = sum(community_sizes) / len(community_sizes)
        print(f'Average community size: {average_size}\n')
    print(f'End comm_tools.community_detection_girvannewman_noweights')
    return Community_list, Num_com_list


def cd_grM_w(G_list):
    Community_list = []
    Num_com_list = []
    print(f'Running comm_tools.cd_greedyModularity_weights')
    for i, g_list in enumerate(G_list):
        # greedy_modularity_communities(G, weight=None, resolution=1, cutoof=1, best_n=None)
        comm_gen = nx.community.greedy_modularity_communities(g_list, weight='weight')  # with weights; runs slower

        communities = [list(community) for community in comm_gen]
        Community_list.append(communities)
        # print(f'Done with iteration:{i}',' \n','Communities detected:', '\n', communities)

    # print(f'-'*120, f'\n Communities:', Community_list)
    for i, c_val in enumerate(Community_list):
        print(f'Communities for jacobian_{i}: ', '\n', c_val, '\n')

        # Compute the number of communities
        num_communities = len(c_val)
        Num_com_list.append(num_communities)
        print(f'Number of communities: {num_communities}')

        # Compute the size of each community in c_val
        community_sizes = [len(community) for community in c_val]
        print(f'Sizes of communities: {community_sizes}')

        # Compute the average size of communities in c_val
        average_size = sum(community_sizes) / len(community_sizes)
        print(f'Average community size: {average_size}\n')
    print(f'End comm_tools.cd_greedyModularity_weights')
    return Community_list, Num_com_list


def cd_grM_nw(G_list):
    Community_list = []
    Num_com_list = []
    print(f'Running comm_tools.cd_greedyModularity_nw')
    for i, g_list in enumerate(G_list):
        # greedy_modularity_communities(G, weight=None, resolution=1, cutoof=1, best_n=None)
        comm_gen = nx.community.greedy_modularity_communities(g_list)  # with weights; runs slower

        communities = [list(community) for community in comm_gen]
        Community_list.append(communities)
        # print(f'Done with iteration:{i}',' \n','Communities detected:', '\n', communities)

    # print(f'-'*120, f'\n Communities:', Community_list)
    for i, c_val in enumerate(Community_list):
        print(f'Communities for jacobian_{i}: ', '\n', c_val, '\n')

        # Compute the number of communities
        num_communities = len(c_val)
        Num_com_list.append(num_communities)
        print(f'Number of communities: {num_communities}')

        # Compute the size of each community in c_val
        community_sizes = [len(community) for community in c_val]
        print(f'Sizes of communities: {community_sizes}')

        # Compute the average size of communities in c_val
        average_size = sum(community_sizes) / len(community_sizes)
        print(f'Average community size: {average_size}\n')
    print(f'End comm_tools.cd_greedyModularity_nw')
    return Community_list, Num_com_list
