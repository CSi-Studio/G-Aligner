import time

import networkx as nx
import numpy as np

from src.map_solver.base_solver import graph_preprocessing, calc_cost_list, flatten_to_matrix_idx, \
    assigned_flatten_idx_to_edges


def greedy_assigner(group_node_list, group_len_list, cost_list):
    assign_time = time.time()
    assigned_idx_list = []
    cost_list[-1] = np.inf
    while np.min(cost_list) != np.inf:
    # for i in range(max(group_len_list)):
        min_cost_idx = np.argmin(cost_list)
        assigned_idx_list.append(min_cost_idx)
        # if i == max(group_len_list) - 1:  # save time for last iteration
        #     break
        idx_tuple = flatten_to_matrix_idx(min_cost_idx, group_len_list)
        settled = np.zeros(group_len_list)
        for j, idx in enumerate(idx_tuple):
            if idx >= len(group_node_list[j]):  # not settle dumb nodes
                continue
            slices = ()
            for k in range(len(idx_tuple)):
                if k == j:
                    slices += idx,
                    continue
                tmp_group_len = group_len_list[k]
                tmp_slice = slice(tmp_group_len)
                slices += tmp_slice,
            settled[slices] = 1
        cost_list[settled.ravel() == 1] = float('inf')
    assign_time = time.time() - assign_time
    return assigned_idx_list, assign_time


def greedy_solve(sub_graph, graph_params, debug=False):
    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]

    # Prepare cost list
    cost_list, cost_time = calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size, greedy=True)

    # Greedy assign
    assigned_idx_list, assign_time = greedy_assigner(group_node_list, group_len_list, cost_list)

    # Confirm assignment result
    assignment_nodes, refine_time = assigned_flatten_idx_to_edges(assigned_idx_list, sub_graph, group_node_list, group_len_list)
    if debug and len(sub_graph) > 20:
        print(len(sub_graph), combination_size, pre_time, cost_time, assign_time, refine_time)
    return assignment_nodes
