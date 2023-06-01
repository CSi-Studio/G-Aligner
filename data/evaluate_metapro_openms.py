import csv
import os

import numpy as np
from pyopenms import *
from src.params import ResultFileReadingParams
from src.result_file_reader import ResultFileReader


def eval(folder_name, align_mz_tolerance, align_rt_tolerance, match_mz_tolerance, match_rt_tolerance):
    tmp_path = os.getcwd()
    result_file_path = os.path.join(tmp_path, folder_name, "metapro")

    result_file_reader = ResultFileReader(
        ResultFileReadingParams(result_file_path, skip_line=0, rt_col_num=2, mz_col_num=1, area_col_num=3))
    result_file_paths, result_file_count = result_file_reader.load_result_paths()

    feature_maps = []
    for path in result_file_paths:
        results = result_file_reader.load_result(path)
        feature_map = FeatureMap()
        for row in results:
            feature = Feature()
            feature.setMZ(row[0])
            feature.setRT(row[1])
            feature.setIntensity(row[2])
            feature_map.push_back(feature)
        feature_maps.append(feature_map)

    # set ref_index to feature map index with the largest number of features
    ref_index = [
        i[0]
        for i in sorted(
            enumerate([fm.size() for fm in feature_maps]), key=lambda x: x[1]
        )
    ][-1]

    aligner = MapAlignmentAlgorithmPoseClustering()
    aligner_params = MapAlignmentAlgorithmPoseClustering().getDefaults()

    aligner_params[b'superimposer:max_shift'] = align_rt_tolerance
    aligner_params[b'superimposer:shift_bucket_size'] = 0.005
    aligner_params[b'superimposer:mz_pair_max_distance'] = align_mz_tolerance
    aligner_params[b'superimposer:rt_pair_distance_fraction'] = 0.005
    aligner_params[b'pairfinder:ignore_charge'] = 'true'
    aligner_params[b'pairfinder:distance_RT:max_difference'] = align_rt_tolerance
    aligner_params[b'pairfinder:distance_MZ:max_difference'] = align_mz_tolerance
    aligner.setReference(feature_maps[ref_index])
    aligner.setParameters(aligner_params)

    # perform alignment and transformation of feature maps to the reference map (exclude reference map)
    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1:]:
        trafo = TransformationDescription()
        aligner.align(feature_map, trafo)
        transformer = MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, True)  # store original RT as meta value

    feature_grouper = FeatureGroupingAlgorithmQT()
    feature_grouper_params = feature_grouper.getDefaults()
    feature_grouper_params[b'ignore_charge'] = 'true'
    feature_grouper_params[b'distance_RT:max_difference'] = match_rt_tolerance
    feature_grouper_params[b'distance_MZ:max_difference'] = match_mz_tolerance
    feature_grouper_params[b'distance_MZ:exponent'] = 1.0
    feature_grouper_params[b'distance_intensity:weight'] = 1.0
    feature_grouper.setParameters(feature_grouper_params)
    consensus_map = ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    # collect information about input maps
    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, ColumnHeader())
        file_description.filename = str(i)
        file_description.size = feature_map.size()
        file_description.unique_id = i
        file_descriptions[i] = file_description

    consensus_map.setColumnHeaders(file_descriptions)
    feature_grouper.group(feature_maps, consensus_map)

    first_line = ['mz', 'rt', 'area', '#']
    for i in range(len(result_file_paths)):
        file_name = os.path.basename(result_file_paths[i]).split('.')[0]
        first_line += [file_name + '_mz', file_name + '_rt', file_name + '_area']

    result_data = np.zeros((consensus_map.size(), 4 + 3 * result_file_count))
    for i in range(consensus_map.size()):
        consensus_feature = consensus_map[i]
        result_data[i, 0] = consensus_feature.getMZ()
        result_data[i, 1] = consensus_feature.getRT()
        result_data[i, 2] = consensus_feature.getIntensity()
        feature_list = consensus_feature.getFeatureList()
        for feature_handle in feature_list:
            idx = feature_handle.getMapIndex()
            mz = feature_handle.getMZ()
            rt = feature_handle.getRT()
            intensity = feature_handle.getIntensity()
            result_data[i, 4 + 3 * idx] = mz
            result_data[i, 5 + 3 * idx] = rt
            result_data[i, 6 + 3 * idx] = intensity
    # ConsensusXMLFile().store('D:\workspace\GAligner\metapro\preview.consensusXML', consensus_map)

    file = open('D:\workspace\GAligner\data\\' + folder_name + '_results_metapro\\' + folder_name + '_aligned_openms.csv', 'w')
    writer = csv.writer(file, dialect='unix', quoting=csv.QUOTE_NONE, quotechar='')
    writer.writerow(first_line)
    writer.writerows(result_data)
    file.close()


eval('TripleTOF_6600', align_mz_tolerance=0.01, align_rt_tolerance=0.5, match_mz_tolerance=0.01, match_rt_tolerance=0.5)
eval('QE_HF', align_mz_tolerance=0.005, align_rt_tolerance=0.3, match_mz_tolerance=0.005, match_rt_tolerance=0.3)
