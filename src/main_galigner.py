
import argparse
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import PolynomialFeatures

from src.params import ResultFileReadingParams, RawFileReadingParams, CoarseRegistrationParams, FineAssignmentParams

from src.raw_file_reader import RawFileReader, Spectrum
from src.result_file_reader import ResultFileReader
from src.coarse_registration import obiwarp_align_3d, ransac_align
from src.fine_alignment import Graph

from src.map_solver.local_bipartite_solver import local_bipartite_solve
from src.map_solver.local_nearest_solver import local_nearest_solve

from src.map_solver.base_solver import graph_partitioning, graph_merging, result_merging
from src.map_solver.gurobi_solver import gurobi_solve
from src.map_solver.ortools_solver import ortools_cp_solve, ortools_mip_solve
from src.map_solver.greedy_solver import greedy_solve
from src.map_solver.vlsns_solver import vlsns_solve

from src.tools import trace_recorder, graph_viewer


class GAligner:
    def __init__(self, result_params, raw_params, coarse_registration_params, fine_assignment_params, saved_folder=None):
        self.result_params = result_params
        self.raw_params = raw_params
        self.coarse_registration_params = coarse_registration_params
        self.fine_assignment_params = fine_assignment_params
        self.debug = False

        if saved_folder is None:
            self.save_folder = str(time.strftime('%Y%m%d_%H%M%S', time.localtime())) + '_' + self.fine_assignment_params.solver
        else:
            self.save_folder = saved_folder
        if self.fine_assignment_params.solver == 'vlsns':
            self.save_folder = self.save_folder + '_' + self.fine_assignment_params.vlsns_solution_init_mode
        trace_recorder.save_params(self.save_folder, result_params, raw_params, coarse_registration_params,
                                   fine_assignment_params)
        print('G-Aligner started...')
        print('Alignment results will be saved in %s' % self.save_folder)

    def plot(self, result):
        plt.rcParams['figure.dpi'] = 400
        for i in range(len(result)):
            idx = (result[i][:, 0] > 15) * (result[i][:, 0] < 20) * (result[i][:, 1] > 502.2) * (result[i][:, 1] < 502.3)
            plt.scatter(result[i][idx, 0], result[i][idx, 1], s=5, marker='.', linewidths=0)
        plt.show()

    def _align(self, result_data_list, result_paths):

        # # 1. Load result data
        # print('Loading result files from disk... ')
        # start_time = time.time()
        # print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time), end='')
        # result_file_reader = ResultFileReader(self.result_params)
        # result_data_list = [result_file_reader.load_result(path) for path in result_paths]
        # print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))
        # print('\tLoaded %d result files.' % len(result_paths))

        result_data_list = [data[np.argsort(data[:, 0])] for data in result_data_list]

        # error = [-0.2, 0.05, -0.15, 0.1, -0.1, 0.15, -0.05, 0.2]
        # for i, result_data in enumerate(result_data_list):
        #     result_data[:, 0] += error[i]

        # 2. Coarse registration
        start_time = time.time()
        print('Performing coarse registration...')
        file_names = []
        for path in result_paths:
            if path is None:
                file_names.append('None')
            else:
                file_names.append(os.path.basename(path))
        print('\t', file_names)
        if self.coarse_registration_params.solver == 'ransac' or (None in result_paths):
            print('\tCalculating warping functions from result data...')
            print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('-' * 50, time.time() - start_time), end='')
            warp_funcs, warp_data = ransac_align(result_data_list, self.coarse_registration_params)
        else:
            print('\tCalculating warping functions from raw data...')
            print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('-' * 50, time.time() - start_time), end='')
            raw_file_reader = RawFileReader(self.raw_params)
            spectra_list = [raw_file_reader.load_ms1_spectra(path) for path in result_paths]
            warp_funcs, warp_data = obiwarp_align_3d(spectra_list, self.coarse_registration_params)
        graph_viewer.plt_coarse_registration_results(warp_data, save=True, save_folder=self.save_folder)
        print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))
        print('\tCoarse aligned figures saved to ~/experiments/%s' % self.save_folder)
        start_time = time.time()
        print('\tApplying warping results to result data...')
        coarse_aligned_result = []  # deep copy
        for result_data in result_data_list:
            coarse_aligned_result.append(result_data.copy())
        for i in range(len(result_data_list)):
            warp_func = warp_funcs[i]
            if warp_func is not None:
                if self.coarse_registration_params.solver == 'ransac' or (None in result_paths):
                    coarse_aligned_result[i][:, 1] = warp_func.predict(
                        PolynomialFeatures(degree=self.coarse_registration_params.degree)
                        .fit_transform(coarse_aligned_result[i][:, 1].reshape(-1, 1)))[:, 0]
                else:
                    coarse_aligned_result[i][:, 1] = warp_func(coarse_aligned_result[i][:, 1])

            done_progress = int((i + 1) / len(result_data_list) * 50)
            print('\r\t[{}{}]{:2d}/{:2d} Time cost:{:.1f}s'
                  .format('▓' * done_progress, '-' * (50 - done_progress), i + 1, len(result_data_list),
                          time.time() - start_time),
                  end='')
        print()

        # 3. Fine assignment
        print('Build graph... ')
        graph = Graph(self.fine_assignment_params)
        sub_graphs, max_nodes = graph.do_build(coarse_aligned_result)
        # sub_graphs, max_nodes = graph.do_build(result_data_list)

        start_time = time.time()
        need_partition = True
        print('Solving MAP for sub-graphs...')
        if self.fine_assignment_params.solver == 'gurobi':
            solver = gurobi_solve
            print('\tUsing gurobi solver...')
        elif self.fine_assignment_params.solver == 'ortools_mip':
            solver = ortools_mip_solve
            print('\tUsing OR-Tools MIP solver...')
        elif self.fine_assignment_params.solver == 'ortools_cp':
            solver = ortools_cp_solve
            print('\tUsing OR-Tools CP solver...')
        elif self.fine_assignment_params.solver == 'greedy':
            solver = greedy_solve
            print('\tUsing greedy solver...')
        elif self.fine_assignment_params.solver == 'vlsns':
            solver = vlsns_solve
            # need_partition = False
            print('\tUsing vlsns solver. init_mode:\'%s\', init_num:%d, update_mode:\'%s\'' %
                  (self.fine_assignment_params.vlsns_solution_init_mode,
                   self.fine_assignment_params.vlsns_solution_init_number,
                   self.fine_assignment_params.vlsns_solution_update_mode))
        elif self.fine_assignment_params.solver == 'local_bipartite':
            solver = local_bipartite_solve
            need_partition = False
            print('\tUsing local bipartite solver...')
        elif self.fine_assignment_params.solver == 'local_nearest':
            solver = local_nearest_solve
            need_partition = False
            print('\tUsing local nearest solver...')
        else:
            return

        subgraph_save_count = 0
        assignment_nodes_list = []
        # get no need to assign number
        from src.map_solver.base_solver import graph_preprocessing
        no_need_to_assign_full_cnt = 0
        no_need_assign_list = []
        need_assign_node_nums = []
        for sub_graph in sub_graphs:
            no_need_to_assign, _, _, _, _ = graph_preprocessing(sub_graph)
            no_need_assign_list.append(no_need_to_assign)
            if no_need_to_assign:
                if len(sub_graph) == len(result_data_list):
                    no_need_to_assign_full_cnt += 1
            else:
                need_assign_node_nums.append(len(sub_graph))
        print('No need to assign number: ' + str(np.sum(no_need_assign_list)))
        print('No need to assign full number: ' + str(no_need_to_assign_full_cnt))
        plt.hist(need_assign_node_nums, bins=list(range(0, max(need_assign_node_nums) + 1, max(1, int(len(result_data_list) / 4)))))
        plt.show()

        need_assign_list = []
        for i, sub_graph in enumerate(sub_graphs):
            if no_need_assign_list[i]:
                assignment_nodes_list.append(sub_graph.nodes(data=True))
                need_assign_list.append(0)
                continue
            node_limit = 20
            if need_partition and len(sub_graph) > node_limit:
                assignment_nodes = []
                sub_sub_graphs = graph_partitioning(sub_graph, node_limit)
                for sub_sub_graph in sub_sub_graphs:
                    sub_assignment_nodes = solver(sub_sub_graph, self.fine_assignment_params, debug=self.debug)
                    assignment_nodes += sub_assignment_nodes
                assignment_nodes = graph_merging(sub_graph, assignment_nodes)
            else:
                assignment_nodes = solver(sub_graph, self.fine_assignment_params, debug=self.debug)
            assignment_nodes_list += assignment_nodes
            need_assign_list += [1] * len(assignment_nodes)

            if len(sub_graph) > len(result_data_list) * 4 and subgraph_save_count < 10:
                subgraph_save_count += 1
                graph_viewer.plt_scatter(sub_graph, save=True, save_folder=self.save_folder, save_name='%d_%d' % (i, len(sub_graph)))
                graph_viewer.plt_assignment(sub_graph, assignment_nodes, save=True, save_folder=self.save_folder,
                                            save_name='%d_%d' % (i, len(sub_graph)))
            done_progress = int((i + 1) / len(sub_graphs) * 50)
            print('\r\t[{}{}] {:d}/{:d} Time cost:{:.1f}s'.format('▓' * done_progress, '-' * (50 - done_progress),
                                                                  i + 1, len(sub_graphs), time.time() - start_time),
                  end='')
        print()
        print('\tSaved figures of %d sub-graphs (nodes > 4 * file_num) to ~/experiments/results/%s'
              % (subgraph_save_count, self.save_folder))

        # assignment_nodes_list = result_merging(assignment_nodes_list,
        #                                        self.fine_assignment_params.mz_tolerance * 2,
        #                                        self.fine_assignment_params.rt_tolerance * 2)
        # 4. Assemble assignment results
        start_time = time.time()
        print('Assembling results...')
        print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('-' * 50, time.time() - start_time), end='')
        result_rows = trace_recorder.prepare_result_rows(result_data_list, assignment_nodes_list, need_assign_list)
        print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))

        return result_rows

    def do_align(self):
        global_start_time = time.time()
        # 1. Load result data and raw data
        print('Finding result files from result folder... ')
        result_file_reader = ResultFileReader(self.result_params)
        result_file_paths, result_file_count = result_file_reader.load_result_paths()
        print('\tFound %d result files.' % result_file_count)

        def recursive_align(file_paths):
            result_data_list = []
            file_names = []
            tmp_result_paths = []
            for path in file_paths:
                if isinstance(path, list):
                    tmp_result = recursive_align(path)
                    result_data_list.append(tmp_result[0])
                    file_names += tmp_result[1]
                    tmp_result_paths.append(None)
                else:
                    file_names.append(os.path.basename(path))
                    result_data_list.append(result_file_reader.load_result(path))
                    tmp_result_paths.append(path)
            return self._align(result_data_list, tmp_result_paths), file_names

        result_data_list, result_file_names = recursive_align(result_file_paths)

        start_time = time.time()
        print('Saving results to file...')
        print('\r\t[{}] 0/ 1 Time cost:{:.1f}s'.format('-' * 50, time.time() - start_time), end='')
        result_file_path = trace_recorder.save_alignment_results(result_data_list, result_file_names, self.save_folder, 4)
        print('\r\t[{}] 1/ 1 Time cost:{:.1f}s'.format('▓' * 50, time.time() - start_time))
        print('\tResult file saved in ' + result_file_path)

        print('G-Aligner execution completed.')
        print('Aligned %d files in %.1fs.' % (len(result_file_names), time.time() - global_start_time))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GAligner parameters')

    # Result file reading params
    parser.add_argument('--result_file_path', type=str, help='Path to result folder', required=True)
    parser.add_argument('--skip_line', type=int, help='Number of header rows ', required=True)
    parser.add_argument('--mz_col_num', type=int, help='M/z column number', required=True)
    parser.add_argument('--rt_col_num', type=int, help='RT column number', required=True)
    parser.add_argument('--area_col_num', type=int, help='Area column number', required=True)

    # Raw file reading params
    parser.add_argument('--raw_file_path', type=str, help='Path to raw file folder', required=True)
    parser.add_argument('--min_intensity', type=float, help='Minimum intensity involved in alignment', required=False, default=10000)

    # Coarse registration parameters (Obiwarp)
    parser.add_argument('--bin_size', type=float, help='', required=False, default=1)
    parser.add_argument('--percent_anchors', type=float, help='', required=False, default=1.0)
    parser.add_argument('--score_type', type=str, help='', required=False, default='pearsons_r2')
    parser.add_argument('--gap_init', type=float, help='', required=False, default=None)
    parser.add_argument('--gap_extend', type=float, help='', required=False, default=None)
    parser.add_argument('--factor_diag', type=float, help='', required=False, default=2)
    parser.add_argument('--factor_gap', type=float, help='', required=False, default=1)
    parser.add_argument('--local_alignment', type=int, help='', required=False, default=0)
    parser.add_argument('--init_penalty', type=float, help='', required=False, default=0)

    # Fine alignment parameters
    parser.add_argument('--mz_tolerance', type=float, help='M/z tolerance', required=True)
    parser.add_argument('--rt_tolerance', type=float, help='RT tolerance', required=True)
    parser.add_argument('--use_ppm', type=bool, help='Use ppm m/z tolerance', required=False, default=False)
    parser.add_argument('--mz_factor', type=float, help='', required=False, default=1)
    parser.add_argument('--rt_factor', type=float, help='', required=False, default=1)
    parser.add_argument('--area_factor', type=float, help='', required=False, default=1)
    # parser.add_argument('', type=, help='', required=False, default=)

    args = parser.parse_args()

    result_file_reading_params = ResultFileReadingParams(result_folder_path=args.result_file_path,
                                                         skip_line=args.skip_line, rt_col_num=args.rt_col_num,
                                                         mz_col_num=args.mz_col_num, area_col_num=args.area_col_num)
    raw_file_reading_params = RawFileReadingParams(raw_file_path=args.raw_file_path, min_intensity=args.min_intensity)

    coarse_registration_params = CoarseRegistrationParams(bin_size=args.bin_size, percent_anchors=args.percent_anchors,
                                                          score_type=args.score_type,
                                                          gap_init=args.gap_init, gap_extend=args.gap_extend,
                                                          factor_diag=args.factor_diag, factor_gap=args.factor_gap,
                                                          local_alignment=args.local_alignment, init_penalty=args.init_penalty)
    fine_assignment_params = FineAssignmentParams(mz_tolerance=args.mz_tolerance, rt_tolerance=args.rt_tolerance,
                                                  use_ppm=args.use_ppm, mz_factor=args.mz_factor,
                                                  rt_factor=args.rt_factor, area_factor=args.area_factor)
    g_aligner = GAligner(result_file_reading_params, raw_file_reading_params, coarse_registration_params,
                         fine_assignment_params)
    g_aligner.do_align()
