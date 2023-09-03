import glob

from pyopenms import ConsensusXMLFile, ConsensusMap, FeatureXMLFile, FeatureMap
import numpy as np
import csv


def convert_consensusxml_to_csv(xml_path, csv_path, file_names):
    consensusMap = ConsensusMap()
    ConsensusXMLFile().load(xml_path, consensusMap)

    csv_file = open(csv_path, 'w')
    writer = csv.writer(csv_file, dialect='unix', quoting=csv.QUOTE_NONE, quotechar='')

    first_row = ['mz', 'rt', 'area', 'need_assign']
    for file_name in file_names:
        first_row += [file_name + "_mz", file_name + "_rt", file_name + "_area"]
    writer.writerow(first_row)

    for i in range(consensusMap.size()):
        consensusFeature = consensusMap[i]
        avg_mz = consensusFeature.getMZ()
        avg_rt = consensusFeature.getRT() / 60
        avg_area = consensusFeature.getIntensity()
        features = np.zeros((len(file_names), 3))
        for feature in consensusFeature.getFeatureList():
            mapIndex = feature.getMapIndex()
            features[mapIndex][0] = feature.getMZ()
            features[mapIndex][1] = feature.getRT() / 60
            features[mapIndex][2] = feature.getIntensity()
        features = features.ravel().tolist()
        row = [avg_mz, avg_rt, avg_area, 0] + features
        writer.writerow(row)
    csv_file.close()
    print(xml_path + " convertion finished.")


def convert_featurexml_to_csv(xml_path, csv_path):
    featureMap = FeatureMap()
    FeatureXMLFile().load(xml_path, featureMap)

    csv_file = open(csv_path, 'w')
    writer = csv.writer(csv_file, dialect='unix', quoting=csv.QUOTE_NONE, quotechar='')
    for feature in featureMap:
        writer.writerow([feature.getMZ(), feature.getRT() / 60, feature.getIntensity()])
    csv_file.close()

    print(xml_path + " convertion finished.")


if __name__ == '__main__':

    wiff_sample_names = ['SampleA_1', 'SampleA_2', 'SampleA_3', 'SampleA_4',
                         'SampleB_1', 'SampleB_2', 'SampleB_3', 'SampleB_4']
    raw_sample_names = ['SA1', 'SA2', 'SA3', 'SA4', 'SA5', 'SB1', 'SB2', 'SB3', 'SB4', 'SB5']

    mtbls_sample_names = ['12W-1', '12W-2', '12W-3', '12W-4', '12W-5', '12W-6', '12W-7', '12W-8',
                          '24W-1', '24W-2', '24W-3', '24W-4', '24W-5', '24W-6', '24W-7', '24W-8',
                          '32W-1', '32W-2', '32W-3', '32W-4', '32W-5', '32W-6', '32W-7', '32W-8',
                          '4W-1', '4W-2', '4W-3', '4W-4', '4W-5', '4W-6', '4W-7', '4W-8',
                          '52W-1', '52W-2', '52W-3', '52W-4', '52W-5', '52W-6', '52W-7', '52W-8']

    # wiff_featurexml_files = glob.glob("D:/workspace/GAligner/data/TripleTOF_6600/openms/*.featureXML")
    # for file in wiff_featurexml_files:
    #     convert_featurexml_to_csv(file, file.replace("featureXML", "csv"))
    #
    # raw_featurexml_files = glob.glob("D:/workspace/GAligner/data/QE_HF/openms/*.featureXML")
    # for file in raw_featurexml_files:
    #     convert_featurexml_to_csv(file, file.replace("featureXML", "csv"))

    mtbls_featurexml_files = glob.glob("D:/workspace/GAligner/data/MTBLS562/openms/*.featureXML")
    for file in mtbls_featurexml_files:
        convert_featurexml_to_csv(file, file.replace("featureXML", "csv"))

    # wiff_consensusxml_file = "D:/workspace/GAligner/data/TripleTOF_6600_results_openms/wiff_aligned.consensusXML"
    # wiff_csv_file = "D:/workspace/GAligner/data/TripleTOF_6600_results_openms/openms_aligned.csv"
    # convert_consensusxml_to_csv(wiff_consensusxml_file, wiff_csv_file, wiff_sample_names)
    #
    # raw_consensusxml_file = "D:/workspace/GAligner/data/QE_HF_results_openms/raw_aligned.consensusXML"
    # raw_csv_file = "D:/workspace/GAligner/data/QE_HF_results_openms/openms_aligned.csv"
    # convert_consensusxml_to_csv(raw_consensusxml_file, raw_csv_file, raw_sample_names)

    mtbls_consensusxml_file = "D:/workspace/GAligner/data/MTBLS562_results_openms/openms_wiff.consensusXML"
    mtbls_csv_file = "D:/workspace/GAligner/data/MTBLS562_results_openms/openms_wiff.csv"
    convert_consensusxml_to_csv(mtbls_consensusxml_file, mtbls_csv_file, mtbls_sample_names)
