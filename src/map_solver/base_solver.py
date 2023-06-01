import time

import numpy as np
import networkx as nx


def graph_preprocessing(graph):
    pre_time = time.time()
    no_need_to_assign = True
    group_node_map = {}
    nodes = list(graph.nodes(data=True))
    areas = []
    for node in nodes:
        areas.append(-node[1]['area'])

    for i in np.argsort(areas):
        node = nodes[i]
        data_idx = node[1]['data_idx']
        if data_idx not in group_node_map.keys():
            group_node_map[data_idx] = [node[0]]
        else:
            no_need_to_assign = False
            group_node_map[data_idx] += [node[0]]
    if no_need_to_assign:
        return no_need_to_assign, None, None, None, None

    group_node_list = list(group_node_map.values())
    group_len_list = []
    for nodes in group_node_list:
        group_len_list.append(len(nodes))
    max_len = max(group_len_list)

    combination_size = 1
    for i in range(len(group_len_list)):
        # add dumb node to each dim
        # if group_len_list[i] < max_len:
        #     group_len_list[i] += 1
        group_len_list[i] += 1
        combination_size *= group_len_list[i]
    pre_time = time.time() - pre_time

    return no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time


def flatten_to_matrix_idx(flatten_idx, group_len_list):
    tmp_result = ()
    for d in range(len(group_len_list) - 1, -1, -1):
        tmp_result += flatten_idx % group_len_list[d],
        flatten_idx = flatten_idx // group_len_list[d]
    return tmp_result[::-1]


def matrix_to_flatten_idx(matrix_idx, group_len_list):
    base = 1
    flatten_idx = 0
    for i in range(len(matrix_idx) - 1, -1, -1):
        if i < len(matrix_idx) - 1:
            base *= group_len_list[i + 1]
        flatten_idx += matrix_idx[i] * base
    return flatten_idx


def flatten_idx_digit_match(flatten_idx, group_len_list, settled_digit, digit_value):
    div = 1
    for group_idx in range(len(group_len_list) - settled_digit, len(group_len_list)):
        div *= group_len_list[group_idx]
    tmp_num = flatten_idx // div
    return tmp_num % group_len_list[len(group_len_list) - 1 - settled_digit] == digit_value


def calc_cost(sub_graph, group_node_list, group_len_list, idx, greedy=False):
    idx_tuple = flatten_to_matrix_idx(idx, group_len_list)
    selected_nodes = []
    dumb_node_num = 0
    for j in range(len(group_len_list)):
        node_idx = idx_tuple[j]
        if node_idx < len(group_node_list[j]):
            selected_nodes.append(group_node_list[j][node_idx])
        else:
            dumb_node_num += 1
    if len(selected_nodes) == 0:
        return -0.7

    # Calculate cost list
    min_rt, max_rt = float('inf'), 0
    min_mz, max_mz = float('inf'), 0
    log_areas = np.zeros(len(selected_nodes))
    for j, node in enumerate(selected_nodes):
        # log_areas[j] = max(0, np.log(sub_graph.nodes[node]['area'] + 1E-6))
        log_areas[j] = sub_graph.nodes[node]['area']
        tmp_mz = sub_graph.nodes[node]['mz']
        tmp_rt = sub_graph.nodes[node]['rt']
        if tmp_mz < min_mz: min_mz = tmp_mz
        if tmp_mz > max_mz: max_mz = tmp_mz
        if tmp_rt < min_rt: min_rt = tmp_rt
        if tmp_rt > max_rt: max_rt = tmp_rt
    # if graph_params.use_ppm:
    #     mz_tolerance = graph_params.mz_tolerance * (min_mz + max_mz) / 2 * 1E-6
    # else:
    #     mz_tolerance = graph_params.mz_tolerance

    # cost
    tmp_sub_graph = sub_graph.subgraph(selected_nodes)
    connected_components = list(nx.connected_components(tmp_sub_graph))

    dumb_node_cost = dumb_node_num * 1.0  # cost of each dumb node should be bigger than the difference of distribution cost of adding a non-dumb new node
    max_adjacent_bonus = 0.3
    # linear_bonus_factor = max_adjacent_bonus / len(group_len_list)
    # max_connected_nodes = max([len(node_set) for node_set in connected_components])
    # node_num_bonus_linear = max_connected_nodes * (max_connected_nodes + 1) / 2 * linear_bonus_factor
    node_num_bonus_linear = 0
    if dumb_node_num == 0 and len(connected_components) == 1:
        node_num_bonus_linear += max_adjacent_bonus

    # if dumb_node_num == len(group_len_list):
    #     node_num_bonus_linear += 0.5

    if greedy:
        node_num_cost = dumb_node_cost - node_num_bonus_linear
    else:
        node_num_cost = -1 * node_num_bonus_linear
    # distribution_cost = graph_params.rt_factor * (max_rt - min_rt) / graph_params.rt_tolerance + \
    #         graph_params.mz_factor * (max_mz - min_mz) / mz_tolerance + \
    #         graph_params.area_factor * (1 - np.min(log_areas) / (np.max(log_areas) + 1E-6))
    #         # graph_params.area_factor * min(1, np.std(log_areas) / (np.mean(log_areas) + 1E-6))
    # distribution_cost /= graph_params.rt_factor + graph_params.mz_factor + graph_params.area_factor
    distribution_cost = sum(nx.get_edge_attributes(nx.minimum_spanning_tree(tmp_sub_graph), 'weight').values())
    # cost = graph_params.rt_factor * np.square((max_rt - min_rt) / graph_params.rt_tolerance) + \
    #         graph_params.mz_factor * np.square((max_mz - min_mz) / mz_tolerance) + \
    #         graph_params.area_factor * np.square(min(1, np.std(log_areas) / (np.mean(log_areas) + 1E-6)))
    # graph_params.area_factor * (1 - min(log_areas)/(max(log_areas) + 1e-6))
    # cost = np.sqrt(cost)
    connectivity_lost = (len(connected_components) - 1) * 0.6

    return node_num_cost + distribution_cost + connectivity_lost


def calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size, greedy=False):
    cost_time = time.time()
    cost_list = np.zeros(combination_size)

    for i in range(combination_size):
        cost_list[i] = calc_cost(sub_graph, group_node_list, group_len_list, i, greedy=greedy)

    cost_time = time.time() - cost_time
    return cost_list, cost_time


def assigned_flatten_idx_to_edges(flatten_idx_list, graph, group_node_list, group_len_list):
    refine_time = time.time()
    assignment_nodes = []
    for flatten_idx in flatten_idx_list:
        tmp_assignment_nodes = []
        idx_tuple = flatten_to_matrix_idx(flatten_idx, group_len_list)
        for j, idx in enumerate(idx_tuple):
            if idx >= len(group_node_list[j]):
                continue
            tmp_assignment_nodes.append(group_node_list[j][idx])
        if len(tmp_assignment_nodes) == 0:
            continue
        # assert len(tmp_assignment_nodes) > 0
        sub_graph = graph.subgraph(tmp_assignment_nodes)
        for node_set in nx.connected_components(sub_graph):
            s = sub_graph.subgraph(node_set)
            assignment_nodes.append(s.nodes(data=True))
    refine_time = time.time() - refine_time
    return assignment_nodes, refine_time


def graph_partitioning(graph):
    def eval(graph, nodes_list):
        cross_edge_num = 0
        for edge in graph.edges:
            for nodes in nodes_list:
                if (edge[0] in nodes and edge[1] not in nodes) or \
                        (edge[1] in nodes and edge[0] not in nodes):
                    cross_edge_num += 1
        return cross_edge_num

    kernighan_lin = nx.algorithms.community.kernighan_lin_bisection(graph)
    # k_clique = nx.algorithms.community.k_clique_communities(graph)
    modularity = nx.algorithms.community.greedy_modularity_communities(graph)
    label_propagation = list(nx.algorithms.community.label_propagation_communities(graph))
    girvan_newman = list(nx.algorithms.community.girvan_newman(graph))
    asyn_fluidc = list(nx.algorithms.community.asyn_fluidc(graph, 2))
    asyn_lpa = list(nx.algorithms.community.asyn_lpa_communities(graph))
    louvain = nx.algorithms.community.louvain_communities(graph)

    kernighan_lin_eval = eval(graph, kernighan_lin)
    modularity_eval = eval(graph, modularity)
    label_propagation_eval = eval(graph, label_propagation)
    girvan_newman_eval = eval(graph, girvan_newman)
    asyn_fluidc_eval = eval(graph, asyn_fluidc)
    asyn_lpa_eval = eval(graph, asyn_lpa)
    louvain_eval = eval(graph, louvain)

    print("kernighan_lin", kernighan_lin_eval, [len(nodes) for nodes in kernighan_lin])
    print("modularity", modularity_eval, [len(nodes) for nodes in modularity])
    print("label_propagation", label_propagation_eval, [len(nodes) for nodes in label_propagation])
    print("girvan_newman", girvan_newman_eval, [len(nodes) for nodes in girvan_newman])
    print("asyn_fluidc", asyn_fluidc_eval, [len(nodes) for nodes in asyn_fluidc])
    print("asyn_lpa", asyn_lpa_eval, [len(nodes) for nodes in asyn_lpa])
    print("louvain", louvain_eval, [len(nodes) for nodes in louvain])
    print("debug")


def graph_partitioning(graph, node_limit):
    sub_graphs = []
    # # communities = nx.algorithms.community.louvain_communities(graph)
    # communities = nx.algorithms.community.greedy_modularity_communities(graph, cutoff=2, resolution=1)
    # # communities = nx.algorithms.community.kernighan_lin_bisection(graph)
    # for nodes in communities:
    #     sub_graph = graph.subgraph(nodes)
    #     if len(nodes) < node_limit:
    #         sub_graphs.append(sub_graph)
    #     else:
    #         sub_graphs += graph_partitioning(sub_graph, node_limit)
    data_idxes = list(nx.get_node_attributes(graph, 'data_idx').values())
    if max(data_idxes) == min(data_idxes):
        return [graph]
    median_idx = np.mean(data_idxes)
    nodes_list = [[], []]
    for node in graph.nodes(data=True):
        if node[1]['data_idx'] < median_idx:
            nodes_list[0].append(node[0])
        else:
            nodes_list[1].append(node[0])
    for nodes in nodes_list:
        sub_graph = graph.subgraph(nodes)
        if len(nodes) <= node_limit:
            sub_graphs.append(sub_graph)
        else:
            sub_graphs += graph_partitioning(sub_graph, node_limit)
    return sub_graphs


def graph_merging(graph, assignment_nodes):
    def data_idx_conflict(data_idxes_1, data_idxes_2):
        for date_idx_1 in data_idxes_1:
            if date_idx_1 in data_idxes_2:
                return True
        return False

    data_idxes = []
    assignment_node_idxes = []
    for nodes in assignment_nodes:
        data_idxes.append([node[1]['data_idx'] for node in nodes])
        assignment_node_idxes.append([node[0] for node in nodes])

    correlations = {}
    for edge in graph.edges:
        edge_idx = [None, None]
        for i, nodes in enumerate(assignment_node_idxes):
            if edge[0] in nodes:
                edge_idx[0] = i
            if edge[1] in nodes:
                edge_idx[1] = i
        assert edge_idx[0] is not None
        assert edge_idx[1] is not None
        if edge_idx[0] == edge_idx[1]:
            continue
        if data_idx_conflict(data_idxes[edge_idx[0]], data_idxes[edge_idx[1]]):
            continue
        edge_idx.sort()
        weight = graph.edges[edge[0], edge[1]]['weight']

        if edge_idx[0] in correlations.keys():
            if edge_idx[1] in correlations[edge_idx[0]].keys():
                if weight < correlations[edge_idx[0]][edge_idx[1]]:
                    correlations[edge_idx[0]][edge_idx[1]] = weight
            else:
                correlations[edge_idx[0]][edge_idx[1]] = weight
        else:
            correlations[edge_idx[0]] = {edge_idx[1]: weight}

    edges = []
    for l in correlations.keys():
        for r in correlations[l].keys():
            edges.append([l, r, correlations[l][r]])

    if len(edges) == 0:
        return assignment_nodes

    sorted_idxes = np.argsort(np.array(edges)[:, -1])

    merge_idxes = []
    for i in sorted_idxes:
        edge = edges[i]
        if data_idx_conflict(data_idxes[edge[0]], data_idxes[edge[1]]):
            continue
        merged_data_idx = data_idxes[edge[0]] + data_idxes[edge[1]]
        matched_idx = [-1, -1]
        for j, idxes in enumerate(merge_idxes):
            if edge[0] in idxes:
                matched_idx[0] = j
            if edge[1] in idxes:
                matched_idx[1] = j

        if matched_idx[0] == -1 and matched_idx[1] == -1:
            merge_idxes.append([edge[0], edge[1]])
            data_idxes[edge[0]] = merged_data_idx
            data_idxes[edge[1]] = merged_data_idx
        elif matched_idx[0] != -1 and matched_idx[1] != -1:
            merge_idxes[matched_idx[0]] += merge_idxes[matched_idx[1]]
            for idx in merge_idxes[matched_idx[0]]:
                data_idxes[idx] = merged_data_idx
            merge_idxes.remove(merge_idxes[matched_idx[1]])
        else:
            merge_idx = matched_idx[0] + matched_idx[1] + 1
            for k in range(2):
                if edge[k] not in merge_idxes[merge_idx]:
                    merge_idxes[merge_idx].append(edge[k])
            for idx in merge_idxes[merge_idx]:
                data_idxes[idx] = merged_data_idx

    remove_idxes = []
    for idxes in merge_idxes:
        idxes.sort()
        assignment_nodes[idxes[0]] = list(assignment_nodes[idxes[0]])
        for j in range(1, len(idxes)):
            assignment_nodes[idxes[0]] += list(assignment_nodes[idxes[j]])
            remove_idxes.append(idxes[j])
    remove_idxes.sort(reverse=True)
    for idx in remove_idxes:
        assignment_nodes.remove(assignment_nodes[idx])

    return assignment_nodes
