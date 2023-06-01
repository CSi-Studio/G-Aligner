import time

import networkx as nx
import gurobipy as gp
from src.map_solver.base_solver import graph_preprocessing, calc_cost_list, flatten_to_matrix_idx,\
    flatten_idx_digit_match, assigned_flatten_idx_to_edges


def gurobi_solve(sub_graph, graph_params, debug=False):
    # Preprocessing
    no_need_to_assign, group_node_list, group_len_list, combination_size, pre_time = graph_preprocessing(sub_graph)
    if no_need_to_assign:
        return [sub_graph.nodes(data=True)]
    # Prepare cost list
    cost_list, cost_time = calc_cost_list(sub_graph, graph_params, group_node_list, group_len_list, combination_size, greedy=True)

    # Set Gurobi optimizer
    assign_time = time.time()
    m = gp.Model("Multidimensional Assignment Problem")
    m.setParam('OutputFlag', 0)

    x = []
    for i in range(combination_size):
        idx_tuple = flatten_to_matrix_idx(i, group_len_list)
        x.append(m.addVar(vtype=gp.GRB.BINARY, name="x" + "%d" * len(group_len_list) % idx_tuple))

    # i: digit from right to left
    for i in range(len(group_len_list)):
        non_dumb_node_num = len(group_node_list[len(group_len_list) - 1 - i])
        for j in range(non_dumb_node_num):
            m.addConstr(gp.quicksum([x[k] for k in range(combination_size)
                                    if flatten_idx_digit_match(k, group_len_list, i, j)]) == 1,
                        name='Constraint%d%d' % (i, j))
    m.addConstr(x[-1] == 0, name='ConstraintDumb')

    m.setObjective(gp.quicksum(x[i] * cost_list[i] for i in range(combination_size)), gp.GRB.MINIMIZE)
    m.optimize()
    assigned_idx_list = []
    for i in range(len(x)):
        if x[i].x > 1e-6:
            assigned_idx_list.append(i)
    assign_time = time.time() - assign_time

    # Confirm assignment result
    assignment_nodes, refine_time = assigned_flatten_idx_to_edges(assigned_idx_list, sub_graph, group_node_list, group_len_list)
    if debug:
        print(len(sub_graph), pre_time, cost_time, assign_time, refine_time)
    return assignment_nodes
