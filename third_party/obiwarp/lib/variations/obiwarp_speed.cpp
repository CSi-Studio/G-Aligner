// STDLIB:
#include <cstdio>
#include <iostream>
#include <fstream>
#include "string.h"

// 3RD PARTY 
#include <argtable2.h>
#include "tnt_stopwatch.h"


// MINE
#include "vec.h"
#include "mat.h"
#include "dynprog.h"
#include "pngio.h"

#define DEBUG (0)


int mymain(const char **SCORE_ARR, int score_cnt);

int main (int argc, char **argv) {
    // Create the argument structures:
    struct arg_str  *score    = arg_str0("s", "score", "<scoretype>", "similarity score to compare vectors");
    struct arg_end  *end      = arg_end(20);

    void* argtable[] = {score,end};
    const char* progname = "obi-warp";
    int nerrors;
    int exitcode=0;

    /* set default values*/
    double diag_factor = 2.f;
    double gap_factor = 1.f;

    /* verify the argtable[] entries were allocated sucessfully */
    if (arg_nullcheck(argtable) != 0)
    {
        /* NULL entries were detected, some allocations must have failed */
        printf("%s: insufficient memory\n",progname);
        exitcode=1;
        goto exit;
    }


    /* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    exitcode = mymain(
            score->sval, score->count
            );

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int mymain(const char **SCORE_ARR, int score_cnt) {
  
    char SCORE[1024];
    if (score_cnt) {
        strcpy(SCORE, SCORE_ARR[0]);
    }
    else {
        strcpy(SCORE, "covariance");
    }

    if (strlen(SCORE) <= 1) {
        strcpy(SCORE, "covariance");
    }


    int num_scores = 10;
    int num_peaks = 1000;

    // Write the file for plotting each guy:
    char plotfilename[1024];
    char plotfilename_toplot[1024];
    strcpy(plotfilename, "speed_test_");
    strcat(plotfilename, SCORE);
    strcpy(plotfilename_toplot, plotfilename);
    strcat(plotfilename_toplot, ".toplot");
    std::ofstream fh(plotfilename_toplot);
    printf("WRITING TO: %s\n", plotfilename_toplot);
    fh << "XYData" << "\n";
    fh << plotfilename << "\n";
    fh << "Scoring Functions Speed Comparison (on draco) " << num_scores << " scores\n";
    fh << "N scans (compared N X N times)\n";
    fh << "sqrt( time to complete " << num_scores << " scores (sec) )\n";
    MatF smat;
    DynProg dyn;

    int i;
    int num_its = 10;
    VecI xaxis(num_its);
    for (i = 0; i < num_its; ++i) {
        xaxis[i] = i * 10;
    }

    fh << SCORE << "numpeaks" << num_peaks << "\n";
    std::cout << "SCORE " << SCORE << "\n";

    VecF yresult(num_its);
    TNT::Stopwatch st;

    for (i = 0; i < xaxis.length(); ++i) {
        int num_scans = xaxis[i];
        MatF mat1(num_scans,num_peaks, 20.f);
        MatF mat2(num_scans,num_peaks, 12.5f);
        MatF smat_slow(mat1.rows(), mat2.rows());

        int cnt1, cnt2;
        VecF *row_vecs1 = new VecF[mat1.rows()];
        VecF *row_vecs2 = new VecF[mat2.rows()];
        mat1.row_vecs(cnt1, row_vecs1);
        mat2.row_vecs(cnt2, row_vecs2);
        if (!strcmp(SCORE, "covariance_slow")) {
            std::cout << "INSIDE" << SCORE << "\n";
            st.start();
            for (int j = 0; j < num_scores; ++j) {
                for (int m = 0; m < cnt1; ++m) {
                    for (int n = 0; n < cnt2; ++n) {
                        smat_slow(m,n) = VecF::covariance(row_vecs1[m], row_vecs2[n]);
                    }
                }
            }
            float timed = st.read();
            yresult[i] = timed;
        }
        else if (!strcmp(SCORE, "pearsonsr_slow")) {
            std::cout << "INSIDE" << SCORE << "\n";
            st.start();
            for (int j = 0; j < num_scores; ++j) {
                for (int m = 0; m < cnt1; ++m) {
                    for (int n = 0; n < cnt2; ++n) {
                        smat_slow(m,n) = VecF::pearsons_r(row_vecs1[m], row_vecs2[n]);
                    }
                }
            }
            float timed = st.read();
            yresult[i] = timed;
        }
        else if (!strcmp(SCORE, "product_slow")) {
            std::cout << "INSIDE" << SCORE << "\n";
            st.start();
            for (int j = 0; j < num_scores; ++j) {
                for (int m = 0; m < cnt1; ++m) {
                    for (int n = 0; n < cnt2; ++n) {
                        smat_slow(m,n) = VecF::dot_product(row_vecs1[m], row_vecs2[n]);
                    }
                }
            }
            float timed = st.read();
            yresult[i] = timed;
        }
        else if (!strcmp(SCORE, "euclidean_slow")) {
            std::cout << "INSIDE" << SCORE << "\n";
            st.start();
            for (int j = 0; j < num_scores; ++j) {
                for (int m = 0; m < cnt1; ++m) {
                    for (int n = 0; n < cnt2; ++n) {
                        smat_slow(m,n) = VecF::euclidean(row_vecs1[m], row_vecs2[n]);
                    }
                }
            }
            float timed = st.read();
            yresult[i] = timed;
        }
        else {
            std::cout << "INSIDE" << SCORE << "\n";
            st.start();
            for (int j = 0; j < num_scores; ++j) {
                dyn.score(mat1, mat2, smat, SCORE, 5);
            }
            float timed = st.read();
            yresult[i] = timed;
        }
    }
    xaxis;
    yresult.square_root();
    xaxis.print(fh,1);
    yresult.print(fh,1);
    fh.close();

    return 0;
}



// Print to file to plot
// title (gap penalty optimization etc.....)
// filename (<b>intercept_linear_gap_penalty_optimization )
    // "slope of gap penalty array"
    // "avg of sq of residuals"
    // mvals ...
    // avgs ...
    // mvals ...
    // avgs ...
    // mvals ...
    // avgs ...

