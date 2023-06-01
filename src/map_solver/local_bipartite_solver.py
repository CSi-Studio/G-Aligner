import time

import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment
from src.map_solver.base_solver import graph_preprocessing


def bipartite_assigner(sub_graph, group_node_list):
    assign_time = time.time()

    node_length_list = [len(nodes) for nodes in group_node_list]
    descend_idxes = np.argsort(-1 * np.array(node_length_list))

    assigned_nodes = group_node_list[descend_idxes[0]].copy()
    assignments = [[node] for node in assigned_nodes]
    for i in range(1, len(descend_idxes)):
        connected_nodes = []
        unconnected_nodes = []
        connected_assignments_idxes = []
        for node in group_node_list[descend_idxes[i]]:
            tmp_assignment_idxes = []
            for j, assignment in enumerate(assignments):
                match_idx = -1
                tmp_min_weight = float('inf')
                for k, assigned_node in enumerate(assignment):
                    if sub_graph.has_edge(node, assigned_node):
                        tmp_weight = sub_graph.get_edge_data(node, assigned_node)['weight']
                        if tmp_weight < tmp_min_weight:
                            tmp_min_weight = tmp_weight
                            match_idx = k
                if match_idx != -1:
                    tmp_assignment_idxes.append([j, match_idx])
            if len(tmp_assignment_idxes) > 0:
                connected_nodes.append(node)
                connected_assignments_idxes.append(tmp_assignment_idxes)
            else:
                unconnected_nodes.append(node)

        if len(connected_nodes) > 0:
            cost_matrix = np.ones([len(connected_nodes), len(assignments)]) * 100
            for m in range(len(connected_nodes)):
                for match_idxes in connected_assignments_idxes[m]:
                    assigned_node = assignments[match_idxes[0]][match_idxes[1]]
                    tmp_weight = sub_graph.get_edge_data(connected_nodes[m], assigned_node)['weight']
                    cost_matrix[m][match_idxes[0]] = tmp_weight
            assigned_idx = linear_sum_assignment(cost_matrix)

            for j in range(len(connected_nodes)):
                if j not in assigned_idx[0]:
                    unconnected_nodes.append(connected_nodes[j])

            for j in range(len(assigned_idx[0])):
                if cost_matrix[assigned_idx[0][j]][assigned_idx[1][j]] == 100:
                    unconnected_nodes.append(connected_nodes[assigned_idx[0][j]])
                    continue
                from_node = connected_nodes[assigned_idx[0][j]]
                assignments[assigned_idx[1][j]].append(from_node)

        for node in unconnected_nodes:
            assignments.append([node])

        assigned_nodes += group_node_list[descend_idxes[i]]

    assign_time = time.time() - assign_time
    return assignments, assign_time


def local_bipartite_solve(sub_graph, graph_params, debug=False):

    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]

    # Local bipartite assign
    assignments, assign_time = bipartite_assigner(sub_graph, group_node_list)

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
