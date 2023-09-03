import time

import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment

from src.map_solver.base_solver import graph_preprocessing, calc_cost, calc_cost_list, flatten_to_matrix_idx, matrix_to_flatten_idx


def init_standard_solution(group_len_list):
    n = max(group_len_list)
    solution = np.zeros([len(group_len_list), n]).astype(np.int32)
    for i in range(len(group_len_list)):
        tmp_len = group_len_list[i]
        solution[i, :tmp_len] = np.array((list(range(tmp_len))))
        solution[i, tmp_len:n] = np.ones(n - tmp_len) * (tmp_len - 1)
    return solution


def init_random_solution(group_len_list):
    n = max(group_len_list)
    solution = np.zeros([len(group_len_list), n]).astype(np.int32)
    for i in range(len(group_len_list)):
        tmp_len = group_len_list[i]
        solution[i, :tmp_len] = np.array((list(range(tmp_len))))
        solution[i, tmp_len:n] = np.ones(n - tmp_len) * (tmp_len - 1)
        np.random.shuffle(solution[i])
    return solution


def init_greedy_solution(greedy_idx_list, group_len_list):
    n = len(greedy_idx_list)
    solution = np.zeros([len(group_len_list), n]).astype(np.int32)
    for i, flatten_idx in enumerate(greedy_idx_list):
        idx_tuple = flatten_to_matrix_idx(flatten_idx, group_len_list)
        solution[:, i] = np.array(list(idx_tuple))
    return solution


def init_local_greedy_solution(assignments, group_node_list):
    solution = np.zeros([len(group_node_list), len(assignments)]).astype(np.int32)
    for i in range(len(group_node_list)):
        solution[i, :] = len(group_node_list[i])
    for i, nodes in enumerate(assignments):
        for node in nodes:
            for j, group_nodes in enumerate(group_node_list):
                if node in group_nodes:
                    solution[j, i] = group_nodes.index(node)
    return solution


def init_random_solutions(group_len_list, solution_count):
    solutions = []
    for i in range(solution_count):
        solutions.append(init_random_solution(group_len_list))
    return solutions


def init_grid_solutions(group_len_list, solution_count):
    solutions = []
    base_solution = init_standard_solution(group_len_list)
    solutions.append(base_solution)
    for i in range(1, solution_count):
        if i % max(group_len_list) == 0:
            base_solution = init_random_solution(group_len_list)
        solution = base_solution.copy()
        for j in range(1, len(group_len_list)):
            solution[j] = np.roll(solution[j], j * i)
        solutions.append(solution)
    return solutions


def vlsns_solve(sub_graph, graph_params, debug=False):
    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]

    # Prepare cost list
    cost_map = {}
    # cost_list, cost_time = calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size)

    # VLSNS assign
    assign_time = time.time()
    init_mode = graph_params.vlsns_solution_init_mode
    init_num = graph_params.vlsns_solution_init_number
    update_mode = graph_params.vlsns_solution_update_mode
    init_solution = None
    if init_mode == 'sss':
        init_solution = init_standard_solution(group_len_list)
    if init_mode == 'ssr':
        init_solution = init_random_solution(group_len_list)
    if init_mode == 'msr':
        init_solution = init_random_solutions(group_len_list, init_num)
    if init_mode == 'msg':
        init_solution = init_grid_solutions(group_len_list, init_num)

    solution, _ = vlsn_search(sub_graph, group_node_list, group_len_list, cost_map, init_solution.copy(), mode=update_mode)
    assign_time = time.time() - assign_time

    def solution_to_assignment(solution, sub_graph, group_node_list):
        assignment_nodes = []
        for i in range(len(solution[0])):
            assignment_idx = solution[:, i]
            tmp_assignment_nodes = []
            for j, idx in enumerate(assignment_idx):
                if idx < len(group_node_list[j]):
                    tmp_assignment_nodes.append(group_node_list[j][idx])
            if len(tmp_assignment_nodes) == 0:
                continue
            tmp_sub_graph = sub_graph.subgraph(tmp_assignment_nodes)
            for node_set in nx.connected_components(tmp_sub_graph):
                s = tmp_sub_graph.subgraph(node_set)
                assignment_nodes.append(s.nodes(data=True))
        return assignment_nodes


    # Confirm assignment result
    refine_time = time.time()

    assignment_nodes = solution_to_assignment(solution, sub_graph, group_node_list)

    # refine_time = time.time() - refine_time
    # if debug:
    #     print(len(sub_graph), pre_time, cost_time, assign_time, refine_time)
    return assignment_nodes


def vlsn_search(sub_graph, group_node_list, group_len_list, cost_map, solution, mode):
    if isinstance(solution, list):  # multi-start mode
        best_solution = None
        minimum_cost = float('inf')
        for tmp_solution in solution:
            optimized_solution, cost = vlsn_search(sub_graph, group_node_list, group_len_list, cost_map, tmp_solution, mode)
            if cost < minimum_cost:
                minimum_cost = cost
                best_solution = optimized_solution
        return best_solution, minimum_cost


    def best_improvement_update(sub_graph, group_node_list, group_len_list, cost_map, solution, taboo_idx):
        improvement_map = {}
        n = len(solution[0])
        original_cost = 0
        for i in range(len(group_len_list)):  # for each dim
            if i == taboo_idx:  # select projection dim
                continue
            cost_matrix = np.zeros([n, n])  # main_row, projection_row
            for m in range(n):
                assignment_idx = solution[:, m].copy()
                for p in range(n):
                    assignment_idx[i] = min(p, group_len_list[i] - 1)
                    cost_idx = matrix_to_flatten_idx(assignment_idx, group_len_list)
                    # cost_matrix[m][p] = cost_list[cost_idx]
                    if cost_idx in cost_map.keys():
                        cost_matrix[m][p] = cost_map[cost_idx]
                    else:
                        cost_matrix[m][p] = calc_cost(sub_graph, group_node_list, group_len_list, cost_idx)
                        cost_map[cost_idx] = cost_matrix[m][p]
            # solve LAP
            original_cost = 0
            for j in range(n):
                original_cost += cost_matrix[j][solution[i][j]]

            optimized_projection_idx = linear_sum_assignment(cost_matrix)[1]

            optimized_cost = 0
            for j in range(n):
                optimized_cost += cost_matrix[j][optimized_projection_idx[j]]
            if optimized_cost < original_cost:
                improvement_map[i] = optimized_cost, optimized_projection_idx
                # solution[i] = optimized_projection_idx
        if len(improvement_map) == 0:
            return None, None, original_cost

        minimum_cost = float('inf')
        update_idx = None
        update_projection_idx = None
        for i in improvement_map.keys():
            if improvement_map[i][0] < minimum_cost:
                minimum_cost = improvement_map[i][0]
                update_idx = i
                update_projection_idx = improvement_map[i][1]
        update_projection_idx[update_projection_idx > group_len_list[update_idx] - 1] = group_len_list[update_idx] - 1
        return update_idx, update_projection_idx, minimum_cost

    def first_improvement_update(sub_graph, group_node_list, group_len_list, cost_map, solution, start_dim_idx):
        n = len(solution[0])
        original_cost = 0
        for i in range(len(group_len_list) - 1):  # for each dim
            dim = (start_dim_idx + i) % len(group_len_list)
            cost_matrix = np.zeros([n, n])  # main_row, projection_row
            for m in range(n):
                assignment_idx = solution[:, m].copy()
                for p in range(n):
                    assignment_idx[dim] = min(p, group_len_list[dim] - 1)
                    cost_idx = matrix_to_flatten_idx(assignment_idx, group_len_list)
                    # cost_matrix[m][p] = cost_list[cost_idx]
                    if cost_idx in cost_map.keys():
                        cost_matrix[m][p] = cost_map[cost_idx]
                    else:
                        cost_matrix[m][p] = calc_cost(sub_graph, group_node_list, group_len_list, cost_idx)
                        cost_map[cost_idx] = cost_matrix[m][p]

            # solve LAP
            original_cost = 0
            for j in range(n):
                original_cost += cost_matrix[j][solution[dim][j]]

            optimized_projection_idx = linear_sum_assignment(cost_matrix)[1]

            optimized_cost = 0
            for j in range(n):
                optimized_cost += cost_matrix[j][optimized_projection_idx[j]]
            if optimized_cost < original_cost:
                optimized_projection_idx[optimized_projection_idx > group_len_list[dim] - 1] = group_len_list[dim] - 1
                return dim, optimized_projection_idx, optimized_cost
        return None, None, original_cost

    taboo_idx = None
    start_dim_idx = 0
    while True:
        if mode == 'best':
            update_idx, update_projection_idx, minimum_cost = best_improvement_update(sub_graph, group_node_list, group_len_list, cost_map, solution, taboo_idx)
            if update_idx is None:
                break
            taboo_idx = update_idx
            solution[update_idx] = update_projection_idx

        if mode == 'first':

            update_idx, update_projection_idx, minimum_cost = first_improvement_update(sub_graph, group_node_list, group_len_list, cost_map, solution, start_dim_idx)
            if update_idx is None:
                break
            solution[update_idx] = update_projection_idx
            start_dim_idx = (update_idx + 1) % len(group_len_list)

    return solution, minimum_cost
