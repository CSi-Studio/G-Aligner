import time

import networkx as nx
import numpy as np

from src.map_solver.base_solver import graph_preprocessing


def nearest_assigner(sub_graph, group_node_list):
    assign_time = time.time()

    node_length_list = [len(nodes) for nodes in group_node_list]
    descend_idxes = np.argsort(-1 * np.array(node_length_list))

    assignments = [[node] for node in group_node_list[descend_idxes[0]]]
    for i in range(1, len(descend_idxes)):
        # build a cost map
        cost_matrix = np.ones([len(group_node_list[descend_idxes[i]]), len(assignments)])
        for j in range(len(group_node_list[descend_idxes[i]])):
            node = group_node_list[descend_idxes[i]][j]
            for k, assignment in enumerate(assignments):
                tmp_min_weight = 100
                for assigned_node in assignment:
                    if sub_graph.has_edge(node, assigned_node):
                        tmp_weight = sub_graph.get_edge_data(node, assigned_node)['weight']
                        if tmp_weight < tmp_min_weight:
                            tmp_min_weight = tmp_weight
                cost_matrix[j][k] = tmp_min_weight

        for j in range(cost_matrix.shape[0]):
            min_idx = np.argmin(cost_matrix)
            min_x_idx = min_idx // cost_matrix.shape[1]
            min_y_idx = min_idx - min_x_idx * cost_matrix.shape[1]
            tmp_node = group_node_list[descend_idxes[i]][min_x_idx]
            tmp_cost = cost_matrix[min_x_idx, min_y_idx]
            if tmp_cost != 100:
                assignments[min_y_idx].append(tmp_node)
            else:
                assignments.append([tmp_node])
            cost_matrix[min_x_idx, :] = float('inf')
            cost_matrix[:, min_y_idx] = float('inf')

    assign_time = time.time() - assign_time
    return assignments, assign_time


def local_nearest_solve(sub_graph, graph_params, debug=False):

    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]

    # Local nearest assign
    assignments, assign_time = nearest_assigner(sub_graph, group_node_list)

    # Confirm assignment result
    start_time = time.time()
    assignment_nodes = []
    for nodes in assignments:
        s = sub_graph.subgraph(nodes)
        assignment_nodes.append(s.nodes(data=True))
    refine_time = time.time() - start_time

    if debug and len(sub_graph) > 20:
        print(len(sub_graph), combination_size, pre_time, assign_time, refine_time)
    return assignment_nodes
