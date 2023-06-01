import os
import time

import networkx as nx
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def plt_coarse_registration_results(coarse_registration_data, save=True, save_folder=''):
    plt.rcParams['savefig.dpi'] = 1200
    plt.figure(figsize=(8, 4))
    path = os.path.join(root_dir, 'experiments', save_folder, 'coarse_alignment_figs')
    if not os.path.exists(path):
        os.mkdir(path)
    # for i in range(len(coarse_registration_data)):
    for i in range(len(coarse_registration_data)):
        ori_rts = coarse_registration_data[i][0]
        warped_rts = coarse_registration_data[i][1]
        plt.plot(warped_rts, warped_rts - ori_rts, label=i+1, linewidth=1)
    plt.xlabel('Warped RT')
    plt.ylabel('Warped RT - Raw RT')
    plt.legend()
    plt.savefig(os.path.join(path, time.strftime('%H%M%S', time.localtime()) + '_coarse_registration_residual.png')) if save else plt.show()
    plt.clf()

    plt.rcParams['savefig.dpi'] = 1200
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    for i in range(len(coarse_registration_data)):
        ori_rts = coarse_registration_data[i][0]
        warped_rts = coarse_registration_data[i][1]
        ints = coarse_registration_data[i][2]
        ax1.plot(ori_rts, ints, label=i+1, linewidth=1)
        ax2.plot(warped_rts, ints, label=i+1, linewidth=1)
    ax1.set_xlabel('Raw RT')
    ax2.set_xlabel('Warped RT')
    ax1.set_ylabel('Intensity')
    ax2.set_ylabel('Intensity')
    plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.45)
    plt.savefig(os.path.join(path, time.strftime('%H%M%S', time.localtime()) + '_coarse_registration_tic.png')) if save else plt.show()
    plt.clf()


def plt_scatter(g, save=False, save_folder='', save_name=''):
    plt.rcParams['savefig.dpi'] = 1200
    mzs = list(nx.get_node_attributes(g, 'mz').values())
    rts = list(nx.get_node_attributes(g, 'rt').values())
    c = list(nx.get_node_attributes(g, 'data_idx'))
    plt.scatter(rts, mzs, c=c, cmap='tab10')
    path = os.path.join(root_dir, 'experiments', save_folder, 'fine_assignment_figs')
    if not os.path.exists(path):
        os.mkdir(path)
    plt.savefig(os.path.join(path, save_name + '_scatter.png')) if save else plt.show()
    plt.clf()


def plt_assignment(g, node_assignment_group, save=False, save_folder='', save_name='', show_weight=True):
    plt.rcParams['savefig.dpi'] = 1200
    plt.rcParams['figure.dpi'] = 400
    mzs = list(nx.get_node_attributes(g, 'mz').values())
    rts = list(nx.get_node_attributes(g, 'rt').values())
    c = list(nx.get_node_attributes(g, 'data_idx'))
    path = os.path.join(root_dir, 'experiments', save_folder, 'fine_assignment_figs')
    if not os.path.exists(path):
        os.mkdir(path)
    plt.scatter(rts, mzs, c=c, s=15, cmap='tab10')
    for assignment_nodes in node_assignment_group:
        color = np.random.rand(3,)
        nodes = []
        for node in assignment_nodes:
            nodes.append(node[0])
        for node_set in nx.connected_components(g.subgraph(nodes)):
            sub_graph = g.subgraph(node_set)
            assignment_edges = nx.minimum_spanning_tree(sub_graph).edges

            for edge in assignment_edges:
                mzs = []
                rts = []
                mzs.append(g.nodes[edge[0]]['mz'])
                mzs.append(g.nodes[edge[1]]['mz'])
                rts.append(g.nodes[edge[0]]['rt'])
                rts.append(g.nodes[edge[1]]['rt'])
                weight = '%.2f' % g.get_edge_data(edge[0], edge[1])['weight']
                plt.plot(rts, mzs, alpha=0.4, c=color)
                if show_weight:
                    plt.text((rts[0] + rts[1]) / 2, (mzs[0] + mzs[1]) / 2, weight)
    plt.savefig(os.path.join(path, save_name + '_assigned.png')) if save else plt.show()
    plt.clf()


def plt_all_edges(g, save=False, save_folder='', save_name='', show_weight=True):
    plt.rcParams['savefig.dpi'] = 400
    plt.rcParams['figure.dpi'] = 400
    mzs = list(nx.get_node_attributes(g, 'mz').values())
    rts = list(nx.get_node_attributes(g, 'rt').values())
    c = list(nx.get_node_attributes(g, 'data_idx'))
    path = os.path.join(root_dir, 'experiments', save_folder, 'fine_assignment_figs')
    if not os.path.exists(path):
        os.mkdir(path)
    plt.scatter(rts, mzs, c=c, s=15, cmap='tab10')

    assignment_edges = g.edges
    for edge in assignment_edges:
        mzs = []
        rts = []
        mzs.append(g.nodes[edge[0]]['mz'])
        mzs.append(g.nodes[edge[1]]['mz'])
        rts.append(g.nodes[edge[0]]['rt'])
        rts.append(g.nodes[edge[1]]['rt'])
        weight = '%.2f' % g.get_edge_data(edge[0], edge[1])['weight']
        plt.plot(rts, mzs, alpha=0.4)
        if show_weight:
            plt.text((rts[0] + rts[1]) / 2, (mzs[0] + mzs[1]) / 2, weight)
    plt.savefig(os.path.join(path, save_name + '_all_edges.png')) if save else plt.show()
    plt.clf()


def show(g):
    pos, labels = _get_layout(g)
    cmap = cm.get_cmap('tab10')
    nx.draw(g, pos=pos, node_size=10, cmap=cmap, node_color=list(labels.values()), font_weight='bold')
    plt.show()


def show_multipartite(g):
    pos, labels = _get_layout(g)
    cmap = cm.get_cmap('tab10')
    pos = nx.multipartite_layout(g, subset_key='data_idx')
    nx.draw(g, pos=pos, node_size=60, cmap=cmap, node_color=list(labels.values()), font_weight='bold')
    plt.show()


def show_assignment(g, assignment_group):
    a = nx.Graph()
    a.add_nodes_from(g.nodes(data=True))
    for assignment in assignment_group:
        a.add_edges_from(assignment)
    show(a)


def _get_layout(g):
    mzs = np.array(list(nx.get_node_attributes(g, 'mz').values()))
    rts = np.array(list(nx.get_node_attributes(g, 'rt').values()))
    pos = {}
    labels = {}
    for node in list(g.nodes(data=True)):
        mz = node[1]['mz']
        rt = node[1]['rt']
        pos_mz = (mz - np.min(mzs))
        pos_rt = (rt - np.min(rts))
        pos[node[0]] = np.array([pos_rt, pos_mz]).astype(np.float32)
        labels[node[0]] = node[1]['data_idx']
    return pos, labels
