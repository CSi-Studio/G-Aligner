import os
import csv
import numpy as np

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def save_params(save_folder, result_file_params, raw_file_params, coarse_registration_params, fine_assignment_params):

    path = os.path.join(root_dir, 'experiments', 'results', save_folder)
    if not os.path.exists(path):
        os.mkdir(path)
    file = open(os.path.join(path, 'params.txt'), 'w')

    # ResultFileReadingParams
    file.write('# Result File Reading Params' + '\n')
    file.write('\tresult_file_path: ' + result_file_params.result_folder_path + '\n')
    file.write('\tskip_line: ' + str(result_file_params.skip_line) + '\n')
    file.write('\trt_col_num: ' + str(result_file_params.rt_col_idx + 1) + '\n')
    file.write('\tmz_col_num: ' + str(result_file_params.mz_col_idx + 1) + '\n')
    file.write('\tarea_col_num: ' + str(result_file_params.area_col_idx + 1) + '\n')
    file.write(os.linesep)

    # RawFileReadingParams
    file.write('# Raw File Reading Params' + '\n')
    file.write('\tmin_intensity: ' + str(raw_file_params.min_intensity) + '\n')
    file.write(os.linesep)

    # CoarseRegistrationParams
    file.write('# Coarse Registration Params' + '\n')
    file.write('\tbin_size: ' + str(coarse_registration_params.bin_size) + '\n')
    file.write('\tpercent_anchors: ' + str(coarse_registration_params.percent_anchors) + '\n')
    file.write('\tscore_type: ' + str(coarse_registration_params.score_type) + '\n')
    file.write('\tgap_init: ' + str(coarse_registration_params.gap_init) + '\n')
    file.write('\tgap_extend: ' + str(coarse_registration_params.gap_extend) + '\n')
    file.write('\tfactor_diag: ' + str(coarse_registration_params.factor_diag) + '\n')
    file.write('\tfactor_gap: ' + str(coarse_registration_params.factor_gap) + '\n')
    file.write('\tlocal_alignment: ' + str(coarse_registration_params.local_alignment) + '\n')
    file.write('\tinit_penalty: ' + str(coarse_registration_params.init_penalty) + '\n')
    file.write(os.linesep)

    # FineAssignmentParams
    file.write('# Fine Assignment Params' + '\n')
    file.write('\trt_tolerance: ' + str(fine_assignment_params.rt_tolerance) + '\n')
    file.write('\tmz_tolerance: ' + str(fine_assignment_params.mz_tolerance) + '\n')
    file.write('\tuse_ppm: ' + str(fine_assignment_params.use_ppm) + '\n')
    file.write('\tmz_factor: ' + str(fine_assignment_params.mz_factor) + '\n')
    file.write('\trt_factor: ' + str(fine_assignment_params.rt_factor) + '\n')
    file.write('\tarea_factor: ' + str(fine_assignment_params.area_factor) + '\n')
    file.write('\tsolver: ' + str(fine_assignment_params.solver) + '\n')
    file.write('\tvlsns_solution_init_mode: ' + str(fine_assignment_params.vlsns_solution_init_mode) + '\n')
    file.write('\tvlsns_solution_init_number: ' + str(fine_assignment_params.vlsns_solution_init_number) + '\n')
    file.write('\tvlsns_solution_update_mode: ' + str(fine_assignment_params.vlsns_solution_update_mode) + '\n')
    file.write(os.linesep)

    file.close()


def prepare_result_rows(result_data_list, assignment_nodes_list, need_assign_list):
    node_start_idxes = [0]
    for i, data in enumerate(result_data_list):
        node_start_idxes += [node_start_idxes[-1] + len(data)]

    row_length = 4
    row_start_idxes = []
    for result_data in result_data_list:
        row_start_idxes.append(row_length)
        tmp_row_length = len(result_data[0])
        if tmp_row_length == 3:
            row_length += 3
        else:
            row_length += tmp_row_length - 4

    rows = []
    for i, nodes in enumerate(assignment_nodes_list):
        assembled_row = np.zeros(row_length)
        tmp_mzs = []
        tmp_rts = []
        tmp_areas = []
        for node in nodes:
            data_idx = node[1]['data_idx']
            row_start_idx = row_start_idxes[data_idx]
            result_data_idx = node[0] - node_start_idxes[data_idx]
            result_data_row = result_data_list[data_idx][result_data_idx]
            if len(result_data_row) == 3:
                assembled_row[row_start_idx: row_start_idx + 3] = np.array(result_data_row)
            else:
                tmp_payload = result_data_row[4:]
                assembled_row[row_start_idx: row_start_idx + len(tmp_payload)] = np.array(tmp_payload)
                need_assign_list[i] = need_assign_list[i] or result_data_row[3]
            tmp_mzs.append(node[1]['mz'])
            tmp_rts.append(node[1]['rt'])
            tmp_areas.append(node[1]['area'])
        tmp_mzs.sort()
        tmp_rts.sort()
        tmp_areas.sort()
        assembled_row[0] = tmp_mzs[len(tmp_mzs) // 2]
        assembled_row[1] = tmp_rts[len(tmp_rts) // 2]
        assembled_row[2] = tmp_areas[len(tmp_areas) // 2]
        assembled_row[3] = need_assign_list[i]
        rows.append(assembled_row)
    return np.array(rows)


def save_alignment_results(result_data_list, file_names, save_folder, min_sample=None):
    if min_sample is None:
        min_sample = (len(file_names) + 1) // 2

    filtered_idxes = np.sum(result_data_list[:, [4 + 3 * i for i in range(len(file_names))]] > 0, axis=-1) > min_sample

    result_root_path = os.path.join(root_dir, 'experiments')
    if not os.path.exists(result_root_path):
        os.mkdir(result_root_path)
    path = os.path.join(result_root_path, save_folder)
    if not os.path.exists(path):
        os.mkdir(path)
    file_path = os.path.join(path, 'aligned_result.csv')
    file = open(file_path, 'w')
    writer = csv.writer(file, dialect='unix', quoting=csv.QUOTE_NONE, quotechar='')
    first_row = ['mz', 'rt', 'area', 'need_assign']
    for file_name in file_names:
        first_row += [file_name + "_mz", file_name + "_rt", file_name + "_area"]
    writer.writerow(first_row)
    writer.writerows(result_data_list[filtered_idxes])
    file.close()
    return file_path
