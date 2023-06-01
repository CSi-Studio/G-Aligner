#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include "string.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "obiwarp/lib/vec.h"
#include "obiwarp/lib/mat.h"
#include "obiwarp/lib/lmat.h"
#include "obiwarp/lib/dynprog.h"


#define DEBUG (0)
namespace py = pybind11;

LMat* create_lmat_from_memory(int len_rt, double *rt, int len_mz, double *mz, double *intensity)
{
    LMat *lmat = new LMat();
    delete lmat->_mz;
    delete lmat->_tm;
    delete lmat->_mat;

    // Get the time values:
    lmat->_tm_vals = len_rt;
    float *tm_tmp = new float[len_rt];
    for(int i=0; i < len_rt; i++) {
      tm_tmp[i] = rt[i];
    }
    lmat->_tm = new VecF(len_rt, tm_tmp);

    // Get the mz values:
    lmat->_mz_vals = len_mz;
    float *mz_tmp = new float[len_mz];
    for(int i=0; i < len_mz; i++) {
      mz_tmp[i] = mz[i];
    }
    lmat->_mz = new VecF(len_mz, mz_tmp);

    // Read the matrix:
    int rows_by_cols = len_rt * len_mz;
    float *mat_tmp = new float[rows_by_cols];
    for(int i=0; i < rows_by_cols; i++) {
      mat_tmp[i] = intensity[i];
    }
    lmat->_mat = new MatF(len_rt, len_mz, mat_tmp);

    return lmat;
}


py::array_t<float> obiwarp(py::array_t<double> py_rt, py::array_t<double> py_mz, py::array_t<double> py_intensity,
                           py::array_t<double> py_rt2, py::array_t<double> py_mz2, py::array_t<double> py_intensity2,
			               float percent_anchors, const char *score,
			               float gap_init, float gap_extend,
			               float factor_diag, float factor_gap,
			               int local_alignment, float init_penalty)
{
    // ************************************************************
    // * CONVERT ARRAY TO MAT
    // ************************************************************
    int len_rt = py_rt.request().size;
    int len_mz = py_mz.request().size;
    int len_rt2 = py_rt2.request().size;
    int len_mz2 = py_mz2.request().size;
    double *rt = (double *)py_rt.request().ptr;
    double *mz = (double *)py_mz.request().ptr;
    double *intensity = (double *)py_intensity.request().ptr;
    double *rt2 = (double *)py_rt2.request().ptr;
    double *mz2 = (double *)py_mz2.request().ptr;
    double *intensity2 = (double *)py_intensity2.request().ptr;
    LMat* lmat1 = create_lmat_from_memory(len_rt, rt, len_mz, mz, intensity);
    LMat* lmat2 = create_lmat_from_memory(len_rt2, rt2, len_mz2, mz2, intensity2);

    // ************************************************************
    // * SCORE THE MATRICES
    // ************************************************************
    if (DEBUG) {
      std::cerr << "Input parameter confirmed!\n";
      std::cerr << " - rt_len = " << lmat1->_tm_vals << "\n";
      std::cerr << " - mz_len = " << lmat1->_mz_vals << "\n";
    }

    MatF smat;
    DynProg dyn;
    dyn.score(*(lmat1->mat()), *(lmat2->mat()), smat, score);

    if (DEBUG) {
      std::cerr << "Matrix scored!\n";
    }

    if (!strcmp(score,"euc")) {
      smat *= -1; // inverting euclidean
    }


    // ************************************************************
    // * PREPARE GAP PENALTY ARRAY
    // ************************************************************

    MatF time_tester;
    MatF time_tester_trans;
    VecF mpt;
    VecF npt;
    VecF mOut_tm;
    VecF nOut_tm;

    int gp_length = smat.rows() + smat.cols();

    VecF gp_array;
    dyn.linear_less_before(gap_extend, gap_init, gp_length, gp_array);

    // ************************************************************
    // * DYNAMIC PROGRAM
    // ************************************************************
    int minimize = 0;
    dyn.find_path(smat, gp_array, minimize, factor_diag, factor_gap, local_alignment, init_penalty);
    if (DEBUG) {
        std::cerr << "Dynamic Time Warping path found!\n";
    }

    VecI mOut;
    VecI nOut;
    dyn.warp_map(mOut, nOut, percent_anchors, minimize);
    if (DEBUG) {
        std::cerr << "Warping anchors decided!\n";
    }

    VecF nOutF;
    VecF mOutF;
    lmat1->tm_axis_vals(mOut, mOutF);
    lmat2->tm_axis_vals(nOut, nOutF);
    lmat2->warp_tm(nOutF, mOutF);
    if (DEBUG) {
        std::cerr << "Piecewise cubic hermite interpolation finished!\n";
    }

    py::array_t<float> warped_rt2 = py::array_t<float>(len_rt2);
    py::buffer_info buffer = warped_rt2.request();

    float* result = (float *)buffer.ptr;
    for(int i=0; i < len_rt2;i++){
      result[i] = lmat2->tm()->pointer()[i];
    }

    delete lmat1;
    delete lmat2;

    return warped_rt2;
}

PYBIND11_MODULE(py_obiwarp, m) {
    m.doc() = "Python Bindings for Obiwarp library.";
    m.attr("__version__") = "0.9.4";

    m.def("obiwarp", &obiwarp, "Perform obiwarp function.");
}

