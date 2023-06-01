import numpy as np
import networkx as nx

from scipy.interpolate import interp1d
from scipy.spatial.distance import cosine

from py_obiwarp import obiwarp
from sklearn.linear_model import RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures

""" OBI-Warp
Python bindings for OBI-Warp algorithm, which was also used in XCMS for mass spectrometry data alignment.

Attributes:
    bin_size: Bin size used in gridding of mass spectrometry data.
    percent_anchors: Percent of anchors used in piecewise cubic hermite interpolation (PCHIP).
    score_type: Score evaluation method used in dynamic time warping (DTW).
                One of ('cor', 'pearsons_r2', 'cov', 'prd', 'euc').
                Default gap penalty of 'mutual_info' was not provided in OBI-Warp paper, not supported here.
    gap_init: Penalty for gap enlargement.
    gap_extend: Penalty for gap enlargement.
    factor_diag: Score weighting factors for diagonal warping in DTW.
    factor_gap: Score weighting factors for gap warping in DTW.
    local_alignment: Perform local alignment instead of the default global alignment.
    init_penalty: Initialized penalty in local alignment.
"""


def _obiwarp(matrix_list, params):
    def score_similarity(rts_1, ints_1, rts_2, ints_2):
        if len(rts_1) != len(ints_1):
            ints_1 = ints_1.reshape((len(rts_1), -1))
            ints_1 = np.sum(ints_1, axis=1)
            ints_2 = ints_2.reshape((len(rts_2), -1))
            ints_2 = np.sum(ints_2, axis=1)

        func_1 = interp1d(rts_1, ints_1, fill_value="extrapolate")
        func_2 = interp1d(rts_2, ints_2, fill_value="extrapolate")
        rt_min = np.floor(min(rts_1[0], rts_2[0]))
        rt_max = np.ceil(max(rts_1[-1], rts_2[-1]))
        min_rt_step = min(np.min(np.diff(rts_1)), np.min(np.diff(rts_2)))
        tick_num = int(np.ceil((rt_max - rt_min) / min_rt_step)) + 1
        rts = np.linspace(rt_min, rt_max, num=tick_num)
        interp_1 = func_1(rts)
        interp_2 = func_2(rts)
        cos_similarity = 1 - cosine(interp_1, interp_2)
        return cos_similarity

    similarity_matrix = np.zeros((len(matrix_list), len(matrix_list)))
    warping_results = []
    for i in range(len(matrix_list)):
        if params.centric_idx != -1 and i != params.centric_idx:
            continue
        warp_func_list = []
        for j in range(len(matrix_list)):
            if i == j:
                warp_func_list.append(None)
                continue
            # init 'percent_anchors'
            # TODO same chrom len ?
            min_rt_len = min(len(matrix_list[i][0]), len(matrix_list[j][0]))
            tmp_percent_anchors = 100 / min_rt_len * 10  # least percent of leaving 10 anchors
            percent_anchors = max(params.percent_anchors, tmp_percent_anchors)
            # percent_anchors = tmp_percent_anchors

            # Perform 1v1 obiwarp (i center)
            if len(matrix_list[i][1]) == 1 and params.score_type != 'euc':
                print('Changing OBI-Warp score_type to \'euc\' for 2D alignment.')
                warped_rts = obiwarp(matrix_list[i][0], matrix_list[i][1], matrix_list[i][2],
                                     matrix_list[j][0], matrix_list[j][1], matrix_list[j][2],
                                     percent_anchors, 'euc', 0.9, 1.8,
                                     params.factor_diag, params.factor_gap, params.local_alignment, params.init_penalty)
            else:
                warped_rts = obiwarp(matrix_list[i][0], matrix_list[i][1], matrix_list[i][2],
                                     matrix_list[j][0], matrix_list[j][1], matrix_list[j][2],
                                     percent_anchors, params.score_type, params.gap_init, params.gap_extend,
                                     params.factor_diag, params.factor_gap, params.local_alignment, params.init_penalty)
            warp_func_list.append(interp1d(matrix_list[j][0], warped_rts, fill_value="extrapolate"))
            similarity = score_similarity(matrix_list[i][0], matrix_list[i][2], warped_rts, matrix_list[j][2])
            similarity_matrix[i][j] = similarity
        warping_results.append(warp_func_list)

    def _mst_warping_results(matrix_list, similarity_matrix, warping_results):
        graph = nx.Graph()
        for i in range(len(similarity_matrix)):
            for j in range(i, len(similarity_matrix)):
                weight = 1 - (similarity_matrix[i][j] + similarity_matrix[j][i]) / 2
                graph.add_edge(i, j, weight=weight)
        mst = nx.minimum_spanning_tree(graph)

        degree_map = np.array(list(nx.degree(mst)))
        center_idx = degree_map[np.argmax(degree_map[:, 1])][0]

        mst_warping_result = []
        for j in range(len(similarity_matrix)):
            if j == center_idx:
                mst_warping_result.append(None)
                continue
            # rts_j = matrix_list[j][0].copy()
            rts = np.linspace(0, matrix_list[j][0][-1], num=int(100*matrix_list[j][0][-1]))
            rts_j = rts.copy()
            node_list = nx.dijkstra_path(mst, center_idx, j)
            for k in range(len(node_list) - 1, 0, -1):
                # Iteratively warp k to the center
                func = warping_results[node_list[k - 1]][node_list[k]]
                rts_j = func(rts_j)
            func_j = interp1d(rts, rts_j, fill_value="extrapolate")
            mst_warping_result.append(func_j)

        return mst_warping_result

    # mst_warping_result = _mst_warping_results(matrix_list, similarity_matrix, warping_results)
    center_idx = np.argmax(np.sum(similarity_matrix, axis=1))
    # center_idx = 0
    mst_warping_result = warping_results[center_idx]
    return mst_warping_result


def obiwarp_align_2d(spectra_list, coarse_registration_params, mode='bpc'):
    def extract_chromatogram(spectra):
        rts = []
        ints = []
        if mode == 'tic':
            for spectrum in spectra:
                rts.append(spectrum.rt)
                ints.append(spectrum.total_ion_current)
        if mode == 'bpc':
            for spectrum in spectra:
                rts.append(spectrum.rt)
                ints.append(spectrum.base_peak_intensity)
        return [np.array(rts).astype(np.float64), [0], np.array(ints).astype(np.float64)]

    chromatograms = []
    for spectra in spectra_list:
        chromatograms.append(extract_chromatogram(spectra))

    return _obiwarp(chromatograms, coarse_registration_params)


def obiwarp_align_3d(spectra_list, coarse_registration_params):
    def _grid_spectra_to_matrix(spectra_list, bin_size):
        # Prepare mzs
        mz_min = float('inf')
        mz_max = 0
        rt_max = 0
        rt_diff_max = 0
        for spectra in spectra_list:
            spectra_rts = []
            for spectrum in spectra:
                spectra_rts.append(spectrum.rt)
                if spectrum.mzs[0] < mz_min:
                    mz_min = np.floor(spectrum.mzs[0])
                if spectrum.mzs[-1] > mz_max:
                    mz_max = np.ceil(spectrum.mzs[-1])
            if spectra[-1].rt > rt_max:
                rt_max = np.ceil(spectra[-1].rt)
            tmp_rt_diff_max = np.max(np.diff(spectra_rts))
            if tmp_rt_diff_max > rt_diff_max:
                rt_diff_max = tmp_rt_diff_max
        mz_tick_num = int(np.ceil((mz_max - mz_min) / bin_size)) + 1
        mz_bins = np.linspace(mz_min, mz_max, num=mz_tick_num)
        # rt_tick_num = int(np.ceil(rt_max / (2 * rt_diff_max)))
        rt_tick_num = 1000
        rt_bins = np.linspace(0, rt_max, num=rt_tick_num)
        # Grid spectra to matrix
        matrix_list = []
        for spectra in spectra_list:
            rts = []
            mzs = []
            ints = []
            for spectrum in spectra:
                rts += [spectrum.rt] * len(spectrum.mzs)
                mzs = np.concatenate((mzs, spectrum.mzs))
                ints = np.concatenate((ints, spectrum.intensities))
                # TODO log
            # find bin max
            Ncount = np.searchsorted(rt_bins, rts, side='right'), np.searchsorted(mz_bins, mzs, side='right')
            rt_on_edge = (rts == rt_bins[-1])
            mz_on_edge = (mzs == mz_bins[-1])
            Ncount[0][rt_on_edge] -= 1
            Ncount[1][mz_on_edge] -= 1
            nbin = np.array([len(rt_bins) + 1, len(mz_bins) + 1])
            xy = np.ravel_multi_index(Ncount, nbin)
            bin_maximum = np.zeros(nbin.prod())
            np.maximum.at(bin_maximum, xy, ints)
            area_max_matrix = bin_maximum.reshape(nbin)
            area_max_matrix = area_max_matrix.astype(float, casting='safe')

            bin_sum = np.bincount(xy, ints, minlength=nbin.prod())
            area_sum_matrix = bin_sum.reshape(nbin)
            area_sum_matrix = area_sum_matrix.astype(float, casting='safe')

            core = (slice(1, -1), slice(1, -1),)
            area_max_matrix = area_max_matrix[core]
            area_sum_matrix = area_sum_matrix[core]

            # int_matrix = np.histogram2d(rts, mzs, [rt_bins, mz_bins], weights=ints, density=True)[0]
            # count_matrix = np.histogram2d(rts, mzs, [rt_bins, mz_bins])[0] + 1e-6
            # mean_matrix = int_matrix / count_matrix
            area_max_matrix = area_max_matrix.ravel(order='C')  # order='C' 'F'
            area_max_matrix = np.log(area_max_matrix + 1E-6)
            matrix_list.append([rt_bins[:-1], mz_bins[:-1], area_max_matrix, area_sum_matrix])
        return matrix_list

    matrix_list = _grid_spectra_to_matrix(spectra_list, coarse_registration_params.bin_size)
    mst_warp_result = _obiwarp(matrix_list, coarse_registration_params)

    coarse_registration_data = []
    for i in range(len(matrix_list)):
        ori_rt = matrix_list[i][0]
        if mst_warp_result[i] is None:
            warped_rt = ori_rt
        else:
            warped_rt = mst_warp_result[i](ori_rt)
        ints = matrix_list[i][3].reshape((len(ori_rt), -1))
        ints = np.max(ints, axis=1)
        coarse_registration_data.append([ori_rt, warped_rt, ints])

    return mst_warp_result, coarse_registration_data


def _ransac(from_data, to_data, params):
    candidate_from_rts = []
    candidate_to_rts = []
    # rt, mz, area
    from_data_idx = 0
    for to_data_idx in range(len(to_data)):
        if to_data[to_data_idx][1] < params.from_rt or to_data[to_data_idx][1] > params.to_rt:
            continue
        if params.use_ppm:
            tmp_mz_tolerance = to_data[to_data_idx][0] * params.mz_tolerance * 1e-6
        else:
            tmp_mz_tolerance = params.mz_tolerance
        tmp_mz_tolerance = tmp_mz_tolerance / 2
        if to_data[to_data_idx][0] < from_data[from_data_idx][0] - tmp_mz_tolerance:
            continue
        while from_data_idx < len(from_data) and from_data[from_data_idx][0] < to_data[to_data_idx][0] - tmp_mz_tolerance:
            from_data_idx += 1
        if from_data_idx >= len(from_data):
            break

        for from_data_iter_idx in range(from_data_idx, len(from_data)):
            if from_data[from_data_iter_idx][0] > to_data[to_data_idx][0] + tmp_mz_tolerance:
                break
            if to_data[to_data_idx][1] - params.rt_tolerance < from_data[from_data_iter_idx][1] < to_data[to_data_idx][1] + params.rt_tolerance:
                candidate_from_rts.append([from_data[from_data_iter_idx][1]])
                candidate_to_rts.append([to_data[to_data_idx][1]])
    candidate_from_rts = np.array(candidate_from_rts)
    candidate_to_rts = np.array(candidate_to_rts)

    # bins = np.linspace(np.min(candidate_from_rts), np.max(candidate_from_rts), 20)
    # bin_ids = np.digitize(candidate_from_rts[:, 0], bins)
    # median_num = max(10, int(np.median([np.sum(bin_ids == i) for i in range(1, 21)])))
    # sampled_ids = []
    # for i in range(1, 21):
    #     np.random.seed(502)
    #     sampled_ids += list(
    #         np.random.choice((bin_ids == i).nonzero()[0], min([median_num, np.sum(bin_ids == i)]), replace=False))

    ransac = RANSACRegressor(residual_threshold=params.rt_residual_threshold, max_trials=1000, random_state=502)
    # ransac.fit(PolynomialFeatures(degree=params.degree).fit_transform(candidate_from_rts[sampled_ids]), candidate_to_rts[sampled_ids])
    ransac.fit(PolynomialFeatures(degree=params.degree).fit_transform(candidate_from_rts), candidate_to_rts)

    # plt_ransac(ransac, candidate_from_rts[sampled_ids], candidate_to_rts[sampled_ids], params.degree)
    # plt_ransac(ransac, candidate_from_rts, candidate_to_rts, params.degree)
    return ransac, ransac.score(PolynomialFeatures(degree=params.degree).fit_transform(candidate_from_rts), candidate_to_rts)


def ransac_align(result_data_list, coarse_registration_params):
    func_map = {}
    score_map = np.zeros((len(result_data_list), len(result_data_list)))
    for i in range(len(result_data_list)):
        if coarse_registration_params.centric_idx != -1 and i != coarse_registration_params.centric_idx:
            continue
        if i not in func_map.keys():
            func_map[i] = {}
        for j in range(len(result_data_list)):
            func, score = _ransac(result_data_list[j], result_data_list[i], coarse_registration_params)
            func_map[i][j] = func
            score_map[i][j] = score

    max_score_idx = np.argmax(np.sum(score_map, axis=-1))
    func_list = []
    warped_data = []
    for i in range(len(result_data_list)):
        ori_rts = result_data_list[i][:, 1]
        sorted_idxes = np.argsort(ori_rts)
        sorted_ori_rts = ori_rts[sorted_idxes]
        sorted_intensities = result_data_list[i][:, 2][sorted_idxes]
        if i == max_score_idx:
            func_list.append(None)
            warped_data.append([sorted_ori_rts, sorted_ori_rts, sorted_intensities])
            continue
        func_list.append(func_map[max_score_idx][i])
        warped_rts = func_list[i].predict(PolynomialFeatures(degree=coarse_registration_params.degree)
                                          .fit_transform(sorted_ori_rts.reshape(-1, 1)))[:, 0]
        warped_data.append([sorted_ori_rts, warped_rts, sorted_intensities])
    return func_list, warped_data


import matplotlib.pyplot as plt


def plt_ransac(func, from_rts, to_rts, degree):
    plt.rcParams['figure.dpi'] = 400
    inlier_mask = func.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    line_x = np.arange(from_rts.min(), from_rts.max()).reshape(-1, 1)
    line_y = func.predict(PolynomialFeatures(degree=degree).fit_transform(line_x))
    # print(line_x[0], line_y[0], line_x[-1], line_y[-1])
    plt.scatter(from_rts[inlier_mask], to_rts[inlier_mask] - from_rts[inlier_mask], s=5, linewidths=0, color='yellowgreen', marker='.', label='Inliers')
    plt.scatter(from_rts[outlier_mask], to_rts[outlier_mask] - from_rts[outlier_mask], s=5, linewidths=0, color='gold', marker='.', label='Outliers')
    plt.plot(line_x[:, 0], line_y[:, 0] - line_x[:, 0], color='cornflowerblue', linewidth=1, label='RANSAC regressor')
    plt.legend(loc='lower right')
    plt.xlabel("Input")
    plt.ylabel("Response")
    plt.show()
