import time

import networkx as nx
import numpy as np

from src.map_solver.base_solver import graph_preprocessing


def local_greedy_assigner(sub_graph):
    assign_time = time.time()

    edge_dist_map = nx.get_edge_attributes(sub_graph, 'weight')
    end_nodes = np.array(list(edge_dist_map.keys()))
    edge_dists = np.array(list(edge_dist_map.values()))
    sorted_idx = np.argsort(edge_dists)

    assignments = []
    assignment_sources = []
    for i in range(len(edge_dists)):
        min_idx = sorted_idx[i]
        nodes = end_nodes[min_idx]
        sources = [sub_graph.nodes[nodes[0]]['data_idx'], sub_graph.nodes[nodes[1]]['data_idx']]

        match_times = 0
        for assignment in assignments:
            for j in range(2):
                if nodes[j] in assignment:
                    match_times += 1
        if match_times == 0:
            assignments.append(list(nodes))
            assignment_sources.append(set(sources))
        if match_times == 1:
            for j in range(len(assignments)):
                if nodes[0] in assignments[j]:
                    if sources[1] not in assignment_sources[j]:
                        assignments[j].append(nodes[1])
                        assignment_sources[j].add(sources[1])
                    else:
                        # make sure no node is disposed
                        assignments.append([nodes[1]])
                        assignment_sources.append({sources[1]})
                    break
                if nodes[1] in assignments[j]:
                    if sources[0] not in assignment_sources[j]:
                        assignments[j].append(nodes[0])
                        assignment_sources[j].add(sources[0])
                    else:
                        assignments.append([nodes[0]])
                        assignment_sources.append({sources[0]})
                    break
        if match_times == 2:
            idxes = []
            for j in range(len(assignments)):
                if nodes[0] in assignments[j]:
                    idxes.append(j)
                if nodes[1] in assignments[j]:
                    idxes.append(j)
            # join assignments
            if idxes[0] != idxes[1] and len(assignment_sources[idxes[0]] & assignment_sources[idxes[1]]) == 0:
                assignments[idxes[0]] += assignments[idxes[1]]
                del assignments[idxes[1]]
                assignment_sources[idxes[0]] = assignment_sources[idxes[0]] | assignment_sources[idxes[1]]
                del assignment_sources[idxes[1]]

    assign_time = time.time() - assign_time
    return assignments, assign_time


def local_greedy_solve(sub_graph, graph_params, debug=False):

    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]

    # Greedy assign
    assignments, assign_time = local_greedy_assigner(sub_graph)

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
