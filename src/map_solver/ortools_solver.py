import time

import networkx as nx
from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model

from src.map_solver.base_solver import graph_preprocessing, calc_cost_list, flatten_to_matrix_idx, \
    flatten_idx_digit_match, assigned_flatten_idx_to_edges


def ortools_mip_solve(sub_graph, graph_params, debug=False):
    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]
    # Prepare cost list
    cost_list, cost_time = calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size)

    assign_time = time.time()

    # Solver
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('PDLP')

    # Variables
    x = []
    for i in range(combination_size):
        idx_tuple = flatten_to_matrix_idx(i, group_len_list)
        x.append(solver.IntVar(0, 1, name="x" + "%d" * len(group_len_list) % idx_tuple))

    # Constraints
    for i in range(len(group_len_list)):
        non_dumb_node_num = len(group_node_list[len(group_len_list) - 1 - i])
        for j in range(non_dumb_node_num):
            solver.Add(solver.Sum([x[k] for k in range(combination_size)
                                   if flatten_idx_digit_match(k, group_len_list, i, j)]) == 1,
                       name='Constraint%d%d' % (i, j))

    # Objective
    objective_terms = []
    for i in range(combination_size):
        objective_terms.append(x[i] * cost_list[i])
    solver.Minimize(solver.Sum(objective_terms))

    # Solve
    status = solver.Solve()

    assigned_idx_list = []
    for i in range(len(x)):
        if x[i].solution_value() > 0.5:
            assigned_idx_list.append(i)
    assign_time = time.time() - assign_time

    # Confirm assignment result
    assignment_nodes, refine_time = assigned_flatten_idx_to_edges(assigned_idx_list, sub_graph,
                                                                                    group_node_list, group_len_list)
    if debug:
        print(len(sub_graph), pre_time, cost_time, assign_time, refine_time)
    return assignment_nodes


def ortools_cp_solve(sub_graph, graph_params, debug=False):
    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]
    # Prepare cost list
    cost_list, cost_time = calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size)

    assign_time = time.time()

    # Model
    model = cp_model.CpModel()

    # Variables
    x = []
    for i in range(combination_size):
        idx_tuple = flatten_to_matrix_idx(i, group_len_list)
        x.append(model.NewBoolVar('x' + '%d' * len(group_len_list) % idx_tuple))

    # Constraints
    for i in range(len(group_len_list)):
        non_dumb_node_num = len(group_node_list[len(group_len_list) - 1 - i])
        for j in range(non_dumb_node_num):
            model.AddExactlyOne([x[k] for k in range(combination_size)
                                 if flatten_idx_digit_match(k, group_len_list, i, j)])

    # Objective
    objective_terms = []
    for i in range(combination_size):
        objective_terms.append(x[i] * cost_list[i])
    model.Minimize(sum(objective_terms))

    # Solve
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    assigned_idx_list = []
    for i in range(len(x)):
        if solver.BooleanValue(x[i]):
            assigned_idx_list.append(i)
    assign_time = time.time() - assign_time

    # Confirm assignment result
    assignment_nodes, refine_time = assigned_flatten_idx_to_edges(assigned_idx_list, sub_graph,
                                                                                    group_node_list, group_len_list)
    if debug:
        print(len(sub_graph), pre_time, cost_time, assign_time, refine_time)
    return assignment_nodes

