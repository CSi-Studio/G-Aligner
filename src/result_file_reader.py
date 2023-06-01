import os
import glob
import time

import numpy as np


class ResultFileReader:
    def __init__(self, result_params):

        if os.path.isfile(result_params.result_folder_path):
            exit('ERROR: result_folder_path must be a folder.')

        self.skip_line = result_params.skip_line
        self.mz_col_idx = result_params.mz_col_idx
        self.rt_col_idx = result_params.rt_col_idx
        self.area_col_idx = result_params.area_col_idx
        self.result_folder_path = result_params.result_folder_path

    def load_result(self, result_path):
        if result_path.endswith('.tsv') or result_path.endswith('.txt') or result_path.endswith('.hills.csv'):
            separator = '\t'
        else:
            separator = ','

        result_file = open(result_path, 'r')
        for j in range(self.skip_line):
            header = result_file.readline().split(separator)
        result_data = np.array([line.strip().replace('"', '').split(separator) for line in result_file])
        results = result_data[:, (self.mz_col_idx, self.rt_col_idx, self.area_col_idx)].astype(np.float32)
        return results

    def _load_result_paths(self, folder_path):
        file_paths = []
        file_count = 0
        files = glob.glob(os.path.join(folder_path, '*'))
        if len(files) == 0:
            return file_paths
        for file in files:
            if os.path.isfile(file) and (file.endswith('.csv') or file.endswith('.tsv') or file.endswith('.txt')):
                file_paths.append(file)
                file_count += 1
            if os.path.isdir(file):
                sub_file_paths, sub_file_count = self._load_result_paths(file)
                if len(sub_file_paths) > 0:
                    file_paths.append(sub_file_paths)
                    file_count += sub_file_count
        return file_paths, file_count

    def load_result_paths(self):
        file_paths, file_count = self._load_result_paths(self.result_folder_path)
        return file_paths, file_count
