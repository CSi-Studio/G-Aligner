
import argparse

from src.params import ResultFileReadingParams, RawFileReadingParams, CoarseRegistrationParams, FineAssignmentParams


def load_params(param_path):
    param_file = open(param_path, 'r')
    lines = param_file.read().split('\n')
    param_map = {}
    for line in lines:
        splits = line.strip().split(':')
        if len(splits) < 2:
            continue
        elif len(splits) == 2:
            param_map[splits[0]] = splits[1]
        elif len(splits) == 3:  # windows path
            param_map[splits[0]] = splits[1] + ':' + splits[2]
        else:
            print('Wrong parameters.')
            return None, None, None, None

    # Result File Reading Params
    result_file_path = param_map['result_file_path']
    skip_line = param_map['skip_line']
    rt_col_num = param_map['rt_col_num']
    mz_col_num = param_map['mz_col_num']
    area_col_num = param_map['area_col_num']

    assert result_file_path is not None
    assert skip_line is not None
    assert rt_col_num is not None
    assert mz_col_num is not None
    assert area_col_num is not None

    result_file_reading_params = ResultFileReadingParams(result_folder_path=result_file_path,
                                                         skip_line=skip_line, rt_col_num=rt_col_num,
                                                         mz_col_num=mz_col_num, area_col_num=area_col_num)

    # Raw File Reading Params
    raw_file_path = param_map['raw_file_path']
    min_intensity = float(param_map['min_intensity'])

    raw_file_reading_params = RawFileReadingParams(raw_file_path=raw_file_path, min_intensity=min_intensity)

    # Coarse Registration Params
    bin_size = float(param_map['bin_size'])
    percent_anchors = float(param_map['percent_anchors'])
    score_type = param_map['score_type']
    gap_init = float(param_map['gap_init'])
    gap_extend = float(param_map['gap_extend'])
    factor_diag = float(param_map['factor_diag'])
    factor_gap = float(param_map['factor_gap'])
    local_alignment = int(param_map['local_alignment'])
    init_penalty = float(param_map['init_penalty'])

    coarse_registration_params = CoarseRegistrationParams(bin_size=bin_size, percent_anchors=percent_anchors,
                                                          score_type=score_type,
                                                          gap_init=gap_init, gap_extend=gap_extend,
                                                          factor_diag=factor_diag, factor_gap=factor_gap,
                                                          local_alignment=local_alignment, init_penalty=init_penalty)

    # Fine Assignment Params
    rt_tolerance = float(param_map['rt_tolerance'])
    mz_tolerance = float(param_map['mz_tolerance'])
    use_ppm = bool(param_map['use_ppm'])
    mz_factor = float(param_map['mz_factor'])
    rt_factor = float(param_map['rt_factor'])
    area_factor = float(param_map['area_factor'])
    solver = param_map['solver']
    vlsns_init_mode = param_map['vlsns_solution_init_mode']
    vlsns_init_num = int(param_map['vlsns_solution_init_number'])
    vlsns_update_mode = param_map['vlsns_solution_update_mode']

    fine_assignment_params = FineAssignmentParams(mz_tolerance=mz_tolerance, rt_tolerance=rt_tolerance,
                                                  use_ppm=use_ppm, mz_factor=mz_factor, rt_factor=rt_factor,
                                                  area_factor=area_factor,
                                                  solver=solver, vlsns_init_mode=vlsns_init_mode,
                                                  vlsns_init_num=vlsns_init_num, vlsns_update_mode=vlsns_update_mode)

    return result_file_reading_params, raw_file_reading_params, coarse_registration_params, fine_assignment_params


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GAligner parameter loader')

    # Result file reading params
    parser.add_argument('--param_path', type=str, help='Path to parameter file', required=True)
    args = parser.parse_args()

    load_params(args.param_path)
    import configparser
    parser = configparser.ConfigParser()
    parser.get()