import glob
import os.path
import warnings

import numpy as np
import pandas as pd

root_path = os.getcwd()

wiff_lib_path = root_path + '/TripleTOF_6600/TripleTOF_6600_annotated.xlsx'
wiff_sample_names = ['SampleA_1', 'SampleA_2', 'SampleA_3', 'SampleA_4',
                     'SampleB_1', 'SampleB_2', 'SampleB_3', 'SampleB_4']

raw_lib_path = root_path + '/QE_HF/QE_HF_annotated.xlsx'
raw_sample_names = ['SA1', 'SA2', 'SA3', 'SA4', 'SA5', 'SB1', 'SB2', 'SB3', 'SB4', 'SB5']

mtbls_lib_path = os.path.join(root_path, 'MTBLS562', 'MTBLS562_annotated.xlsx')
mtbls_sample_names = ['12W-1', '12W-2', '12W-3', '12W-4', '12W-5', '12W-6', '12W-7', '12W-8',
                      '24W-1', '24W-2', '24W-3', '24W-4', '24W-5', '24W-6', '24W-7', '24W-8',
                      '32W-1', '32W-2', '32W-3', '32W-4', '32W-5', '32W-6', '32W-7', '32W-8',
                      '4W-1', '4W-2', '4W-3', '4W-4', '4W-5', '4W-6', '4W-7', '4W-8',
                      '52W-1', '52W-2', '52W-3', '52W-4', '52W-5', '52W-6', '52W-7', '52W-8']

raw_result_path = []


def load_lib(lib_path, sample_names):
    warnings.simplefilter('ignore')
    data = pd.read_excel(lib_path)
    ana_idx = data['Type'] == 'ANALYTES'
    lib_matrix = []
    for i, sample in enumerate(sample_names):
        mzs = np.array(data[sample + '_Mz'][ana_idx])
        rts = np.array(data[sample + '_Rt'][ana_idx])
        lib_matrix.append(np.vstack([mzs, rts]).transpose())
    lib_matrix = np.array(lib_matrix)
    sorted_idx = np.argsort(np.mean(lib_matrix[:, :, 0], axis=0))
    lib_matrix = lib_matrix[:, sorted_idx, :]
    return lib_matrix

def load_lib_dev(lib_path, sample_names):
    warnings.simplefilter('ignore')
    data = pd.read_excel(lib_path)
    ana_idx = data['Type'] == 'ANALYTES'
    lib_matrix = []
    lib_mzs = np.array(data['m/z'][ana_idx])
    lib_rts = np.array(data['RT'][ana_idx])
    for i, sample in enumerate(sample_names):
        mzs = np.array(data[sample + '_Mz'][ana_idx])
        rts = np.array(data[sample + '_Rt'][ana_idx])
        mz_devs = np.array(data[sample + '_Mz'][ana_idx]) - lib_mzs
        rt_devs = np.array(data[sample + '_Rt'][ana_idx]) - lib_rts
        lib_matrix.append(np.vstack([mzs, rts, mz_devs, rt_devs]).transpose())
    lib_matrix = np.array(lib_matrix)
    sorted_idx = np.argsort(lib_mzs)
    lib_matrix = lib_matrix[:, sorted_idx, :]
    return lib_matrix


def load_results(result_paths, separator, skip_line, rt_col_idx, mz_col_idx, area_col_idx):
    result_matrix = []
    for result_path in result_paths:
        result_file = open(result_path, 'r')
        for i in range(skip_line):
            header = result_file.readline().split(separator)
        result_data = np.array([line.strip().split(separator) for line in result_file])
        mzs = result_data[:, mz_col_idx].astype(np.float32)
        rts = result_data[:, rt_col_idx].astype(np.float32)
        areas = result_data[:, area_col_idx].astype(np.float32)
        sorted_idx = np.argsort(mzs)
        result = np.vstack([mzs[sorted_idx], rts[sorted_idx], areas[sorted_idx]]).transpose()
        available_idx = [True] * len(result)
        for i in range(1, len(result)):
            if result[i][0] == result[i - 1][0] and result[i][1] == result[i - 1][1]:
                # find max
                available_idx[i] = False
        result_matrix.append(result[available_idx])

    return result_matrix


def load_aligned_result(result_path, error=None, separator=',', skip_line=1):
    result_file = open(result_path, 'r')
    for i in range(skip_line):
        header = result_file.readline().split(separator)
    aligned_matrix = np.array([line.strip().split(separator) for line in result_file]).astype(np.float32)
    aligned_matrix = aligned_matrix[np.argsort(aligned_matrix[:, 0])]
    if error is not None:
        for i in range(len(error)):
            aligned_matrix[:, 5 + 3 * i] -= error[i]
    return aligned_matrix


# def load_aligned_result(result_path, area_col_idxes, separator=',', skip_line=1, rt_col_idx=1, mz_col_idx=0):
#     result_file = open(result_path, 'r')
#     for i in range(skip_line):
#         header = result_file.readline().split(separator)
#     result_data = np.array([line.strip().split(separator) for line in result_file])
#     rts = result_data[:, rt_col_idx].astype(np.float32)
#     mzs = result_data[:, mz_col_idx].astype(np.float32)
#     areas = result_data[:, area_col_idxes].astype(np.float32)
#     sorted_idx = np.argsort(mzs)
#     aligned_matrix = np.concatenate([rts[sorted_idx].reshape(-1, 1), mzs[sorted_idx].reshape(-1, 1), areas[sorted_idx]], axis=1)
#     return aligned_matrix


def match_result_to_lib(lib_matrix, result_matrix, mz_tolerance, use_ppm, rt_tolerance):
    match_matrix = np.ones(lib_matrix.shape[:2]).astype(np.int32) * -1
    mz_error_matrix = np.zeros(lib_matrix.shape[:2])
    rt_error_matrix = np.zeros(lib_matrix.shape[:2])
    for i in range(len(lib_matrix)):
        lib_data = lib_matrix[i]
        sorted_lib_idx = np.argsort(lib_data[:, 0])
        result_data = result_matrix[i]
        result_idx = 0
        for j in range(len(lib_data)):
            lib_idx = sorted_lib_idx[j]
            tmp_mz_tolerance = mz_tolerance
            if use_ppm:
                tmp_mz_tolerance = lib_data[lib_idx][0] * mz_tolerance * 1e-6
            mz_start = lib_data[lib_idx][0] - tmp_mz_tolerance
            mz_end = lib_data[lib_idx][0] + tmp_mz_tolerance
            if result_data[result_idx][0] > mz_end:
                continue
            while result_idx < len(result_data) and result_data[result_idx][0] < mz_start:
                result_idx += 1
            if result_idx >= len(result_data):
                break

            rt_start = lib_data[lib_idx][1] - rt_tolerance
            rt_end = lib_data[lib_idx][1] + rt_tolerance
            match_idx = -1
            min_dist = float('inf')
            mz_error = 0
            rt_error = 0
            for result_iter_idx in range(result_idx, len(result_data)):
                if result_data[result_iter_idx][0] > mz_end:
                    break
                if (result_data[result_iter_idx][1] < rt_start) \
                        or (result_data[result_iter_idx][1] > rt_end):
                    continue
                tmp_rt_error = abs(result_data[result_iter_idx][1] - lib_data[lib_idx][1])
                tmp_mz_error = abs(result_data[result_iter_idx][0] - lib_data[lib_idx][0])
                dist = np.square(tmp_rt_error / rt_tolerance) + np.square(tmp_mz_error / mz_tolerance)
                if dist < min_dist:
                    min_dist = dist
                    mz_error = tmp_mz_error
                    rt_error = tmp_rt_error
                    match_idx = result_iter_idx
            match_matrix[i, lib_idx] = match_idx
            mz_error_matrix[i, lib_idx] = mz_error
            rt_error_matrix[i, lib_idx] = rt_error
    return match_matrix, mz_error_matrix, rt_error_matrix


def match_aligned_to_result(aligned_matrix, result_matrix):
    match_matrix = np.ones((aligned_matrix.shape[0], int((aligned_matrix.shape[1] - 4) / 3))).astype(np.int32) * -1
    for i in range(len(result_matrix)):
        sample_result = result_matrix[i]
        aligned_result = aligned_matrix[:, [4 + 3 * i, 5 + 3 * i, 6 + 3 * i]]
        result_idx = 0
        for aligned_idx in range(len(aligned_result)):
            aligned_mz = aligned_result[aligned_idx][0]
            aligned_rt = aligned_result[aligned_idx][1]
            aligned_area = aligned_result[aligned_idx][2]
            if aligned_area < 1e-6:
                continue
            while sample_result[result_idx][0] < aligned_mz - 0.1:
                result_idx += 1
            for result_iter_idx in range(result_idx, len(sample_result)):
                if abs(aligned_mz - sample_result[result_iter_idx][0]) < 1e-6 and \
                        abs(aligned_rt - sample_result[result_iter_idx][1]) < 1e-6 and \
                        abs(aligned_area - sample_result[result_iter_idx][2]) < 1e-6:
                    match_matrix[aligned_idx, i] = result_iter_idx
                    break
                if sample_result[result_iter_idx][0] > aligned_mz + 0.1:
                    break
    return match_matrix


def match_aligned_area_to_result(aligned_matrix, result_matrix, mz_tolerance, use_ppm, rt_tolerance):
    match_matrix = np.ones((aligned_matrix.shape[0], int((aligned_matrix.shape[1] - 4) / 3))).astype(np.int32) * -1
    for i in range(len(result_matrix)):
        sample_result = result_matrix[i]
        aligned_result = aligned_matrix[:, [4 + 3 * i, 5 + 3 * i, 6 + 3 * i]]
        sorted_aligned_idx = np.argsort(aligned_result[:, 0])
        result_idx = 0
        for j in range(len(aligned_result)):
            aligned_idx = sorted_aligned_idx[j]
            target_area = aligned_result[aligned_idx][2]
            if target_area < 1e-6:
                continue
            tmp_mz_tolerance = mz_tolerance
            if use_ppm:
                tmp_mz_tolerance = aligned_result[aligned_idx][0] * mz_tolerance * 1e-6
            mz_start = aligned_result[aligned_idx][0] - tmp_mz_tolerance
            mz_end = aligned_result[aligned_idx][0] + tmp_mz_tolerance
            if sample_result[result_idx][0] > mz_end:
                continue
            while result_idx < len(sample_result) and sample_result[result_idx][0] < mz_start:
                result_idx += 1
            if result_idx >= len(sample_result):
                break

            rt_start = aligned_result[aligned_idx][1] - rt_tolerance
            rt_end = aligned_result[aligned_idx][1] + rt_tolerance
            for result_iter_idx in range(result_idx, len(sample_result)):
                if sample_result[result_iter_idx][0] > mz_end:
                    break
                if (sample_result[result_iter_idx][1] < rt_start) \
                        or (sample_result[result_iter_idx][1] > rt_end):
                    continue
                if abs(sample_result[result_iter_idx][2] - target_area) < 1e-6:
                    match_matrix[aligned_idx, i] = result_iter_idx
                    break
    return match_matrix


def eval_alignment_performance(lib_match_matrix, align_match_matrix, need_assign_list):
    eval_result = []
    for i in range(len(lib_match_matrix[0])):
        # 970
        lib_match_assignment = lib_match_matrix[:, i]
        if sum(lib_match_assignment != -1) < len(lib_match_matrix) / 2:
            lib_match_assignment = np.ones(len(lib_match_matrix)) * -1
        matched_idxes = []
        for j in range(len(lib_match_matrix)):
            if lib_match_assignment[j] == -1:
                continue
            matched_idx = np.where(align_match_matrix[:, j] == lib_match_assignment[j])[0]
            if len(matched_idx) > 0:
                matched_idxes.append(matched_idx[0])
        if len(matched_idxes) > 0:
            max_idx = np.argmax(np.bincount(matched_idxes))
            max_assignment = align_match_matrix[max_idx]
            need_assign = need_assign_list[max_idx]
        else:
            max_assignment = np.ones(len(lib_match_matrix)).astype(np.int32) * -1
            need_assign = False
        tp = np.sum((max_assignment == lib_match_assignment) * (max_assignment != -1))
        fp = np.sum(max_assignment[max_assignment != -1] != lib_match_assignment[max_assignment != -1])
        tn = np.sum((max_assignment == -1) * (lib_match_assignment == -1))
        fn = np.sum((max_assignment == -1) * (lib_match_assignment != -1))
        tt = np.sum((max_assignment != -1) * (lib_match_assignment == -1))
        eval_result.append([tp, fp, tn, fn, tt, need_assign])
    eval_result = np.array(eval_result)
    def eval_performance(eval_result):
        total_accuracy = np.sum(eval_result[:, [0, 2]]) / np.sum(eval_result[:, 0:4])
        total_precision = np.sum(eval_result[:, 0]) / np.sum(eval_result[:, [0, 1]])
        total_recall = np.sum(eval_result[:, 0]) / np.sum(eval_result[:, [0, 3]])
        total_f1 = 2 * total_precision * total_recall / (total_precision + total_recall)
        comp_acc = np.sum(eval_result[:, 0] + eval_result[:, 2] == len(lib_match_matrix)) / len(eval_result)
        total_scores = [total_accuracy, total_precision, total_recall, total_f1, comp_acc]
        print(np.sum(eval_result, axis=0)[:4].tolist() + total_scores)
    eval_performance(eval_result)
    # eval_performance(eval_result[eval_result[:, -1] == 1])
    return eval_result


def eval_galigner(lib_matrix, result_matrix, aligned_matrix, mz_tolerance=0.0001, rt_tolerance=0.0001):

    align_area_idxes = [6 + 3 * i for i in range(len(lib_matrix))]
    need_assign_list = aligned_matrix[:, 3] == 1
    aligned_area = aligned_matrix[:, align_area_idxes]
    result_num = np.sum([len(result) for result in result_matrix])
    feature_num = np.sum(aligned_area != 0)
    miss_num = len(aligned_matrix) * len(lib_matrix) - feature_num
    # print("Results:", len(aligned_matrix), result_num, feature_num, miss_num,
    #       feature_num / result_num, feature_num / (feature_num + miss_num))
    # print(np.mean(np.sum(aligned_area != 0, axis=-1)),
    #       np.sum(np.sum(aligned_area != 0, axis=-1) == len(lib_matrix)),
    #       np.sum(np.sum(aligned_area != 0, axis=-1) == len(lib_matrix)) / len(aligned_matrix))
    lib_match_matrix, mz_error_matrix, rt_error_matrix = match_result_to_lib(lib_matrix, result_matrix, mz_tolerance=mz_tolerance, use_ppm=False, rt_tolerance=rt_tolerance)
    align_match_matrix = match_aligned_area_to_result(aligned_matrix, result_matrix, mz_tolerance=0.0001, use_ppm=False, rt_tolerance=0.5)
    eval_result = eval_alignment_performance(lib_match_matrix, align_match_matrix, need_assign_list)
    # print()

    # match_number_list = np.sum(lib_match_matrix >= 0, axis=0)
    # plt.rcParams['figure.dpi'] = 400
    # n, bins, patches = plt.hist(match_number_list, bins=len(result_paths), range=(1, len(lib_sample_names)+1), rwidth=0.8, align='left')
    # for i in range(len(n)):
    #     plt.text(bins[i], n[i] + max(n) * 0.01, int(n[i]), fontsize=10, horizontalalignment="center")
    # plt.show()
    # plt.scatter(mz_error_matrix[0][lib_match_matrix[0] >= 0], rt_error_matrix[0][lib_match_matrix[0] >= 0], s=3)
    # plt.show()
    # # plt.scatter(mz_error_matrix[0][match_idx_matrix[0] >= 0] / lib_matrix[0][match_idx_matrix[0] >= 0][:, 1] * 1e6,
    # #             rt_error_matrix[0][match_idx_matrix[0] >= 0], s=3)
    # # plt.show()
    # plt.hist2d(mz_error_matrix[0][lib_match_matrix[0] >= 0], rt_error_matrix[0][lib_match_matrix[0] >= 0])
    # plt.show()


if __name__ == '__main__':

    result_root_path = os.getcwd()

    wiff_result_paths = [glob.glob(os.path.join(result_root_path, 'TripleTOF_6600', 'metapro',
                                                '*' + name + '*'))[0] for name in wiff_sample_names]
    raw_result_paths = [glob.glob(os.path.join(result_root_path, 'QE_HF', 'metapro',
                                               '*' + name + '*'))[0] for name in raw_sample_names]
    mtbls_result_paths = [glob.glob(os.path.join(result_root_path, 'MTBLS562', 'metapro',
                                               name + '*'))[0] for name in mtbls_sample_names]

    wiff_lib_matrix = load_lib(wiff_lib_path, wiff_sample_names)
    raw_lib_matrix = load_lib(raw_lib_path, raw_sample_names)
    mtbls_lib_matrix = load_lib(mtbls_lib_path, mtbls_sample_names)
    mtbls_lib_matrix_dev = load_lib_dev(mtbls_lib_path, mtbls_sample_names)

    wiff_result_matrix = load_results(wiff_result_paths, separator=',', skip_line=0, rt_col_idx=1, mz_col_idx=0, area_col_idx=2)
    raw_result_matrix = load_results(raw_result_paths, separator=',', skip_line=0, rt_col_idx=1, mz_col_idx=0, area_col_idx=2)
    mtbls_result_matrix = load_results(mtbls_result_paths, separator=',', skip_line=0, rt_col_idx=1, mz_col_idx=0, area_col_idx=2)

    # print('wiff')
    # aligned_paths = [os.path.join(path, 'aligned_result.csv') for path in glob.glob(os.path.join(result_root_path, 'TripleTOF_6600_results_metapro', '[0-9]*'))]
    # aligned_paths += [os.path.join(result_root_path, 'TripleTOF_6600_results_metapro', 'TripleTOF_6600_aligned_mzmine2_ransac.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'TripleTOF_6600_results_metapro', 'TripleTOF_6600_aligned_openms.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'TripleTOF_6600_results_metapro', 'TripleTOF_6600_group_aligned_xcms.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'TripleTOF_6600_results_metapro', 'TripleTOF_6600_obiwarp_aligned_xcms.csv')]
    #
    # for aligned_result_path in aligned_paths:
    #     print(aligned_result_path)
    #     aligned_matrix = load_aligned_result(aligned_result_path)
    #     eval_galigner(wiff_lib_matrix, wiff_result_matrix, aligned_matrix)
    #
    #
    # print('raw')
    # aligned_paths = [os.path.join(path, 'aligned_result.csv') for path in glob.glob(os.path.join(result_root_path,
    #                                                                                              'QE_HF_results_metapro', '[0-9]*'))]
    # aligned_paths += [os.path.join(result_root_path, 'QE_HF_results_metapro', 'QE_HF_aligned_mzmine2_ransac.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'QE_HF_results_metapro', 'QE_HF_aligned_openms.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'QE_HF_results_metapro', 'QE_HF_group_aligned_xcms.csv')]
    # aligned_paths += [os.path.join(result_root_path, 'QE_HF_results_metapro', 'QE_HF_obiwarp_aligned_xcms.csv')]
    #
    # for aligned_result_path in aligned_paths:
    #     print(aligned_result_path)
    #     aligned_matrix = load_aligned_result(aligned_result_path)
    #     eval_galigner(raw_lib_matrix, raw_result_matrix, aligned_matrix)

    print('mtbls')
    aligned_paths = [os.path.join(path, 'aligned_result.csv') for path in glob.glob(os.path.join(result_root_path,
                                                                                                 'MTBLS562_results_metapro', '[0-9]*'))]
    aligned_paths += [os.path.join(result_root_path, 'MTBLS562_results_metapro', 'MTBLS562_aligned_mzmine2_join.csv')]
    aligned_paths += [os.path.join(result_root_path, 'MTBLS562_results_metapro', 'MTBLS562_aligned_mzmine2_ransac.csv')]
    aligned_paths += [os.path.join(result_root_path, 'MTBLS562_results_metapro', 'MTBLS562_aligned_openms.csv')]
    aligned_paths += [os.path.join(result_root_path, 'MTBLS562_results_metapro', 'MTBLS562_group_aligned_xcms.csv')]
    aligned_paths += [os.path.join(result_root_path, 'MTBLS562_results_metapro', 'MTBLS562_obiwarp_aligned_xcms.csv')]

    for aligned_result_path in aligned_paths:
        print(aligned_result_path)
        aligned_matrix = load_aligned_result(aligned_result_path)
        eval_galigner(mtbls_lib_matrix, mtbls_result_matrix, aligned_matrix)

