import sys
import time
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


class Graph:
    def __init__(self, graph_params):
        self.mz_tolerance = graph_params.mz_tolerance
        self.rt_tolerance = graph_params.rt_tolerance
        self.use_ppm = graph_params.use_ppm
        self.mz_factor = graph_params.mz_factor
        self.rt_factor = graph_params.rt_factor
        self.area_factor = graph_params.area_factor

    def do_build(self, data_list):
        g = nx.Graph()
        data_lens = []
        start_idxes = [0]
        node_idx = 0

        print('\tConverting result to graph nodes...')
        start_time = time.time()
        print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('-' * 50, time.time() - start_time), end='')
        for i, data in enumerate(data_list):
            data_lens += [len(data)]
            start_idxes += [start_idxes[-1] + len(data)]
            for point in data:
                g.add_node(node_idx, data_idx=i, mz=point[0], rt=point[1], area=point[2])
                node_idx += 1
        print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))

        print('\tAnalyzing node relations and build edges...')
        start_time = time.time()
        comparison_idx = 0
        total_comparisons = ((len(data_list) - 1) // 2) * len(data_list) + ((len(data_list) - 1) % 2) * (
                len(data_list) // 2)

        def build_edges(g, left_data, right_data, left_start_idx, right_start_idx, mz_tolerance, use_ppm, rt_tolerance,
                        mz_factor, rt_factor, area_factor, window_factor):
            edges = []
            right_idx = 0
            for left_idx in range(len(left_data)):
                adj_nodes = np.array(list(nx.edges(g, left_start_idx + left_idx)))
                if len(adj_nodes) > 0 and np.sum((adj_nodes[:, 1] >= right_start_idx) *
                                                 (adj_nodes[:, 1] < right_start_idx + len(right_data))) > 0:
                    continue
                tmp_mz_tolerance = mz_tolerance
                if use_ppm:
                    tmp_mz_tolerance = left_data[left_idx][0] * mz_tolerance * 1e-6
                mz_start = left_data[left_idx][0] - tmp_mz_tolerance * window_factor
                mz_end = left_data[left_idx][0] + tmp_mz_tolerance * window_factor
                if right_data[right_idx][0] > mz_end:
                    continue
                while right_idx < len(right_data) and right_data[right_idx][0] < mz_start:
                    right_idx += 1
                if right_idx >= len(right_data):
                    break

                rt_start = left_data[left_idx][1] - rt_tolerance * window_factor
                rt_end = left_data[left_idx][1] + rt_tolerance * window_factor
                for right_iter_idx in range(right_idx, len(right_data)):
                    if right_data[right_iter_idx][0] > mz_end:
                        break
                    if (right_data[right_iter_idx][1] < rt_start) \
                            or (right_data[right_iter_idx][1] > rt_end):
                        continue
                    # dist = mz_factor * abs(left_data[left_idx][0] - right_data[right_iter_idx][0]) / tmp_mz_tolerance + \
                    #     rt_factor * abs(left_data[left_idx][1] - right_data[right_iter_idx][1]) / rt_tolerance + \
                    #     area_factor * (1 - min(left_data[left_idx][2], right_data[right_iter_idx][2]) /
                    #                    (max(left_data[left_idx][2], right_data[right_iter_idx][2]) + 1E-6))
                    # dist /= mz_factor + rt_factor + area_factor
                    dist = mz_factor * abs(left_data[left_idx][0] - right_data[right_iter_idx][0]) / tmp_mz_tolerance + \
                        rt_factor * abs(left_data[left_idx][1] - right_data[right_iter_idx][1]) / rt_tolerance
                    dist /= mz_factor + rt_factor
                    edges.append([left_start_idx + left_idx, right_start_idx + right_iter_idx, dist])
            edges = np.array(edges)
            # filtered_edges = np.array([False] * len(edges))
            # for left_idx in set(edges[:, 0]):
            #     edge_idx = edges[:, 0] == left_idx
            #     min_dist = np.min(edges[edge_idx][:, 2])
            #     filtered_edges[edge_idx * ((edges[:, 2] < 2 * min_dist) + (edges[:, 2] < min_dist + 0.2))] = True
            # for right_idx in set(edges[:, 1]):
            #     edge_idx = edges[:, 1] == right_idx
            #     min_dist = np.min(edges[edge_idx][:, 2])
            #     filtered_edges[edge_idx * ((edges[:, 2] < 2 * min_dist) + (edges[:, 2] < min_dist + 0.2))] = True

            # for edge in edges[filtered_edges]:
            for edge in edges:
                g.add_edge(int(edge[0]), int(edge[1]), weight=edge[2])

        for i in range(len(data_list)):
            for j in range(i + 1, len(data_list)):
                left_data = data_list[i]
                right_data = data_list[j]
                # build_edges(g, left_data, right_data, start_idxes[i], start_idxes[j], self.mz_tolerance, self.use_ppm,
                #             self.rt_tolerance / 2, self.mz_factor, self.rt_factor, self.area_factor, pow(2, -3))
                # build_edges(g, left_data, right_data, start_idxes[i], start_idxes[j], self.mz_tolerance, self.use_ppm,
                #             self.rt_tolerance, self.mz_factor, self.rt_factor, self.area_factor, pow(2, -2))
                # build_edges(g, left_data, right_data, start_idxes[i], start_idxes[j], self.mz_tolerance, self.use_ppm,
                #             self.rt_tolerance / 2, self.mz_factor, self.rt_factor, self.area_factor, pow(2, -1))
                build_edges(g, left_data, right_data, start_idxes[i], start_idxes[j], self.mz_tolerance, self.use_ppm,
                            self.rt_tolerance, self.mz_factor, self.rt_factor, self.area_factor, 1)
                comparison_idx += 1
                done_progress = int((comparison_idx / total_comparisons) * 50)
                print('\r\t[{}{}]{:2d}/{:2d} Time cost:{:.1f}s'.format('▓' * done_progress, '-' * (50 - done_progress),
                                                                       comparison_idx, total_comparisons,
                                                                       time.time() - start_time), end='')
        print()
        print('\tSplitting graph into sub-graphs...')
        start_time = time.time()
        print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time), end='')
        sub_graphs = [g.subgraph(c) for c in nx.connected_components(g)]
        max_nodes = 0
        max_edges = 0
        for s in sub_graphs:
            if max_nodes < len(s.nodes):
                max_nodes = len(s.nodes)
            if max_edges < len(s.edges):
                max_edges = len(s.edges)
        print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))
        print('\tGraph built. Total %d sub-graphs. Max %d nodes, %d edges.' % (len(sub_graphs), max_nodes, max_edges))

        for sub_graph in sub_graphs:
            max_area_map = {}
            nodes = sub_graph.nodes(data=True)
            for node in nodes:
                if (node[1]['data_idx'] not in max_area_map.keys()) or (max_area_map[node[1]['data_idx']] < node[1]['area']):
                    max_area_map[node[1]['data_idx']] = node[1]['area']
            for edge in sub_graph.edges:
                dist = sub_graph.edges[edge]['weight'] * (self.mz_factor + self.rt_factor)
                area_0 = nodes[edge[0]]['area'] / (max_area_map[nodes[edge[0]]['data_idx']] + 1e-6)
                area_1 = nodes[edge[1]]['area'] / (max_area_map[nodes[edge[1]]['data_idx']] + 1e-6)
                dist += self.area_factor * (1 - min(area_0, area_1) / (max(area_0, area_1) + 1e-6))
                dist /= self.mz_factor + self.rt_factor + self.area_factor
                sub_graph.edges[edge]['weight'] = dist
        return sub_graphs, max_nodes
