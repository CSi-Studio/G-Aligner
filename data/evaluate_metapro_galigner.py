import os
from src.params import ResultFileReadingParams, RawFileReadingParams, CoarseRegistrationParams, FineAssignmentParams
from src.main_galigner import GAligner


def eval_wiff(folder_name, solver, vlsns_init_mode):
    result_file_path = os.path.join(os.getcwd(), folder_name, 'metapro')

    result_file_reading_params = ResultFileReadingParams(result_file_path, skip_line=0,
                                                         rt_col_num=2, mz_col_num=1, area_col_num=3)
    raw_file_reading_params = RawFileReadingParams()
    coarse_registration_params = CoarseRegistrationParams(solver='ransac', mz_tolerance=0.01, use_ppm=False, centric_idx=0,
                                                          rt_tolerance=0.5, rt_residual_threshold=0.05, degree=1)
    fine_assignment_params = FineAssignmentParams(rt_tolerance=0.3, mz_tolerance=0.01, mz_factor=1, rt_factor=1,
                                                  area_factor=1, use_ppm=False, solver=solver)
    if vlsns_init_mode is not None:
        fine_assignment_params.vlsns_solution_init_mode = vlsns_init_mode

    g_aligner = GAligner(result_file_reading_params, raw_file_reading_params, coarse_registration_params,
                         fine_assignment_params)
    g_aligner.do_align()


def eval_raw(folder_name, solver, vlsns_init_mode):
    result_file_path = os.path.join(os.getcwd(), folder_name, 'metapro')

    result_file_reading_params = ResultFileReadingParams(result_file_path, skip_line=0,
                                                         rt_col_num=2, mz_col_num=1, area_col_num=3)
    raw_file_reading_params = RawFileReadingParams()
    coarse_registration_params = CoarseRegistrationParams(solver='ransac', mz_tolerance=0.005, use_ppm=False, centric_idx=0,
                                                          rt_tolerance=0.3, rt_residual_threshold=0.02, degree=1)
    fine_assignment_params = FineAssignmentParams(rt_tolerance=0.15, mz_tolerance=0.003, mz_factor=1, rt_factor=1,
                                                  area_factor=1, use_ppm=False, solver=solver)
    if vlsns_init_mode is not None:
        fine_assignment_params.vlsns_solution_init_mode = vlsns_init_mode

    g_aligner = GAligner(result_file_reading_params, raw_file_reading_params, coarse_registration_params,
                         fine_assignment_params)
    g_aligner.do_align()


# TripleTOF_6600
eval_wiff('TripleTOF_6600', solver='local_bipartite', vlsns_init_mode=None)
eval_wiff('TripleTOF_6600', solver='greedy', vlsns_init_mode=None)
eval_wiff('TripleTOF_6600', solver='gurobi', vlsns_init_mode=None)
eval_wiff('TripleTOF_6600', solver='vlsns', vlsns_init_mode='msr')
eval_wiff('TripleTOF_6600', solver='vlsns', vlsns_init_mode='msg')


# QE_HF
eval_raw('QE_HF', solver='local_bipartite', vlsns_init_mode=None)
eval_raw('QE_HF', solver='greedy', vlsns_init_mode=None)
eval_raw('QE_HF', solver='gurobi', vlsns_init_mode=None)
eval_raw('QE_HF', solver='vlsns', vlsns_init_mode='msr')
eval_raw('QE_HF', solver='vlsns', vlsns_init_mode='msg')
