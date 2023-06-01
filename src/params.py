
class ResultFileReadingParams:
    def __init__(self, result_folder_path, skip_line, rt_col_num, mz_col_num, area_col_num):

        self.result_folder_path = result_folder_path
        self.skip_line = skip_line
        self.rt_col_idx = rt_col_num - 1
        self.mz_col_idx = mz_col_num - 1
        self.area_col_idx = area_col_num - 1


class RawFileReadingParams:
    def __init__(self, min_intensity=10000):
        self.min_intensity = min_intensity


class CoarseRegistrationParams:
    def __init__(self, solver='ransac', centric_idx=-1, from_rt=0, to_rt=float('inf'), mz_tolerance=5, use_ppm=True, rt_tolerance=0.5, rt_residual_threshold=0.05, degree=1,
                 bin_size=1, percent_anchors=1.0, score_type='pearsons_r2', gap_init=None, gap_extend=None,
                 factor_diag=2, factor_gap=1, local_alignment=0, init_penalty=0):

        self.solver = solver
        self.centric_idx = centric_idx

        # RANSAC params
        self.from_rt = from_rt
        self.to_rt = to_rt
        self.mz_tolerance = mz_tolerance
        self.use_ppm = use_ppm
        self.rt_tolerance = rt_tolerance
        self.rt_residual_threshold = rt_residual_threshold
        self.degree = degree

        # Obi-Warp params
        self.bin_size = bin_size
        self.percent_anchors = percent_anchors
        self.score_type = score_type
        self.gap_init = gap_init
        self.gap_extend = gap_extend
        self.factor_diag = factor_diag
        self.factor_gap = factor_gap
        self.local_alignment = local_alignment
        self.init_penalty = init_penalty

        # Parameter completion
        if gap_init is None:
            if self.score_type in ['cor', 'pearsons_r2']: self.gap_init = 0.3
            if self.score_type == 'cov': self.gap_init = 0.0
            if self.score_type == 'prd': self.gap_init = 0.0
            if self.score_type == 'euc': self.gap_init = 0.9
        if gap_extend is None:
            if self.score_type in ['cor', 'pearsons_r2']: self.gap_extend = 2.4
            if self.score_type == 'cov': self.gap_extend = 11.7
            if self.score_type == 'prd': self.gap_extend = 7.8
            if self.score_type == 'euc': self.gap_extend = 1.8

        # Parameter verification
        assert solver in ['ransac', 'obiwarp']
        assert isinstance(self.mz_tolerance, float) or isinstance(self.mz_tolerance, int)
        assert isinstance(self.use_ppm, bool)
        assert isinstance(self.rt_tolerance, float) or isinstance(self.rt_tolerance, int)
        assert isinstance(self.rt_residual_threshold, float) or isinstance(self.rt_residual_threshold, int)

        assert isinstance(self.bin_size, float) or isinstance(self.bin_size, int)
        assert isinstance(percent_anchors, float) or isinstance(percent_anchors, int) and 0 <= percent_anchors <= 100
        assert score_type in ['cor', 'pearsons_r2', 'cov', 'prd', 'euc']
        assert isinstance(self.gap_init, float) or isinstance(self.gap_init, int)
        assert isinstance(self.gap_extend, float) or isinstance(self.gap_extend, int)
        assert isinstance(self.factor_diag, float) or isinstance(self.factor_diag, int)
        assert isinstance(self.factor_gap, float) or isinstance(self.factor_gap, int)
        assert isinstance(self.local_alignment, int)
        assert isinstance(self.init_penalty, float) or isinstance(self.init_penalty, int)


class FineAssignmentParams:
    def __init__(self, rt_tolerance=0.2, mz_tolerance=0.005, use_ppm=False, mz_factor=1, rt_factor=1, area_factor=1, solver='vlsns',
                 vlsns_init_mode='ssg', vlsns_init_num=20, vlsns_update_mode='best'):
        self.rt_tolerance = rt_tolerance
        self.mz_tolerance = mz_tolerance
        self.use_ppm = use_ppm
        self.mz_factor = mz_factor
        self.rt_factor = rt_factor
        self.area_factor = area_factor

        """ Fine alignment solvers
        MAP (multi-dimensional assignment problem) solvers
            'gurobi': An all-purpose standard linear solver, commercial, need license
            'ortools_mip': Google OR-Tools solver for Mixed Integer Programming (MIP)
            'ortools_cp': Google OR-Tools solver for Cone Programming (CP)
            'greedy': An optimization method using assignment cost greedy algorithm
            'vlsns': An optimization method using VLSNS (Very Large-Scale Neighborhood Search)
        Local optimum solvers
            'local_greedy': An optimization method using node distance greedy algorithm
            'local_bipartite': An optimization method based on bipartite assignment algorithm
            'local_nearest': An optimization method based on nearest search
        """
        self.solver = solver

        """ VLSNS Params
            Solution initialization modes
                'sss': Single-Start with Standard solution
                'ssr': Single-Start with Random solution
                'ssg': Single-Start with Greedy solution
                'msr': Multi-Start with Random solutions
                'msg': Multi-Start with Gridded solutions
                
            Solution initialization number
                number of initial solutions, used in 'msr' and 'msg'
                
            Solution update modes
                'best': Update the dimension, which has the BEST cost descent in optimization iteration.
                'first': Update the dimension, which is the FIRST cost descent in optimization iteration.
        """
        self.vlsns_solution_init_mode = vlsns_init_mode
        self.vlsns_solution_init_number = vlsns_init_num
        self.vlsns_solution_update_mode = vlsns_update_mode

        # Parameter verification
        assert isinstance(rt_tolerance, float) or isinstance(rt_tolerance, int)
        assert isinstance(mz_tolerance, float) or isinstance(mz_tolerance, int)
        assert isinstance(use_ppm, bool)
        assert isinstance(rt_factor, float) or isinstance(rt_factor, int) and rt_factor >= 0
        assert isinstance(mz_factor, float) or isinstance(mz_factor, int) and mz_factor >= 0
        assert isinstance(area_factor, float) or isinstance(area_factor, int) and area_factor >= 0
        assert mz_factor + rt_factor + area_factor > 0, 'Wrong factors. Sum of factors cannot be 0.'

        assert solver in ['gurobi', 'ortools_mip', 'ortools_cp', 'greedy', 'vlsns', 'local_greedy', 'local_bipartite',
                          'local_nearest'],\
            'Wrong fine alignment solver.' \
            ' Select within [\'gurobi\', \'ortools_mip\', \'ortools_cp\', \'greedy\', \'vlsns\',' \
            ' \'local_greedy\', \'local_bipartite\', \'local_nearest\']'
        assert vlsns_init_mode in ['sss', 'ssr', 'ssg', 'msr', 'msg'], \
            'Wrong VLSNS init mode. Select within [\'sss\', \'ssr\', \'ssg\', \'msr\', \'msg\'].\n' \
            '\t\'sss\': Single-Start with Standard solution.\n' \
            '\t\'ssr\': Single-Start with Random solution.\n' \
            '\t\'ssg\': Single-Start with Greedy solution.\n' \
            '\t\'msr\': Multi-Start with Random solutions.\n' \
            '\t\'msg\': Multi-Start with Gridded solutions.'
        assert isinstance(vlsns_init_num, int) and vlsns_init_num > 0, \
            'Illegal VLSNS init number. Input a number bigger than 0.'
        assert vlsns_update_mode in ['best', 'first'], \
            'Wrong VLSNS update mode. Select within [\'best\', \'first\'].\n' \
            '\t\'best\': Update the dimension, which has the BEST cost descent in optimization iteration.\n' \
            '\t\'first\': Update the dimension, which is the FIRST cost descent in optimization iteration.'
