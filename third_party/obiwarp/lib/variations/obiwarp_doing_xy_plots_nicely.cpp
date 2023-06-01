// STDLIB:
#include <cstdio>
#include <iostream>
#include <fstream>
#include "string.h"

// 3RD PARTY 
#include <argtable2.h>


// MINE
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "dynprog.h"
#include "pngio.h"

#define DEBUG (0)


int mymain(int AXES, int BINARY, const char **logs, int nlogs, int IMAGES, int LOCAL,
        const char **outfiles, int outfile_cnt, const char **SCORE_ARR, int score_cnt,
        double diag_factor, double gap_factor,
        const char **testfile, int testfile_cnt,
        const char **infiles, int ninfiles,
        const char **smat_out_files, int nsmat_out_files,
        const char **smat_in_files, int nsmat_in_files);

int main (int argc, char **argv) {
    // Create the argument structures:
    struct arg_lit  *axes     = arg_lit0("a", "axes", "file contains x and y axis labels");
    struct arg_lit  *binary   = arg_lit0("b", "binary", "binary data file");
    struct arg_str  *logbase  = arg_str0(NULL, "log", "<base>", "takes the log of the smat");
    struct arg_lit  *help     = arg_lit0("h", "help", "prints this help and exits");
    struct arg_lit  *images   = arg_lit0("i", "images", "png images created");
    struct arg_lit  *local    = arg_lit0("l", "local", "local alignment");
    struct arg_file *outfiles = arg_file0("o", "outfile", "<outfile>", "output to file");
    struct arg_str  *score    = arg_str0("s", "score", "<scoretype>", "similarity score to compare vectors");
    struct arg_dbl  *diagfac  = arg_dbl0("d", "diagfac", "2.0", "factor applied to diagonal steps");
    struct arg_dbl  *gapfac   = arg_dbl0("g", "gapfac", "1.0", "factor applied to off-diagonal steps");
    struct arg_file *testfile = arg_file0("t", "testfile", "<times.txt>", "file with times to verify alignment");
    struct arg_file *infiles  = arg_filen(NULL, NULL, NULL, 2, 2, "files to align (first is template)");
    struct arg_file *smat_out = arg_file0(NULL, "smat_out", "<file>", "writes binary smat to file and exits");
    struct arg_file *smat_in  = arg_file0(NULL, "smat_in", "<file>", "uses binary smat file instead of calculating");
    struct arg_end  *end      = arg_end(20);

    void* argtable[] = {axes,binary,logbase,help,images,local,outfiles,score,diagfac,gapfac,testfile,infiles,smat_out,smat_in,end};
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

    if (help->count > 0 || argc < 3)
    {
//        "*****************************************************************\n" <<
//        "usage: obi-warp [-a] [-b] [-i] [-l] [-t <olap.txt>] [-s <scoretype>] file1 file2\n" <<
//        "*****************************************************************\n" <<
        printf("******************************************************************************\n");
        printf("Usage: %s ", progname);
        arg_print_syntax(stdout,argtable,"\n");
        printf("\n");
        printf("Interpolated dynamic time warping\n");
        printf("\n");
        arg_print_glossary(stdout,argtable,"  %-26s %s\n");
        printf("\n");
        printf("\n");
        printf("FILE FORMAT: \n");
        printf("Data should be in an m(rows)x n(cols) matrix (space delimited)\n");
        printf("where each line contains one row of data.  Should be same # cols.\n"); 
        printf("Data will be aligned along the m axis.\n");
        printf("ARGUMENTS (default marked by asterik*): \n");
        printf("\n");
        printf("-b, --binary = file is binary rather than *ascii\n");
        printf("     binary format: all vals 4 bytes\n"); 
        printf("     (int) # rows, [IF axes:] (floats) m axis values,\n");
        printf("     (int) # cols, [IF axes:] (floats) n axis values,\n");
        printf("     (floats) matrix data values (packed together) row1, row2, row3 ...\n");
        printf("\n");
        printf("-a, --axes = 1st line in file contains num x coordinates, 2nd the coords\n"); 
        printf("     3rd line in file contains num y coordinates, 4th, the coords\n"); 
        printf("     [space delimited]\n");
        printf("\n");
        printf("-s, --score = scoring function: *covariance, product (dot product)\n"); 
        printf("     pearsons_r, pearsons_r2, euclidean, mutual_info\n");
        printf("\n");
        printf("-l, --local = local rather than *global alignment\n");
        printf("\n");
        printf("-i, --images = creates png images of the alignment process\n");
        printf("     (if png2matrix installed)\n");
        printf("\n");
        printf("-t, --test = given a file of two columns of time points,\n");
        printf("     emit the avg sum of sqrs of residuals and avg time diff\n");
        printf("******************************************************************************\n");

        exitcode=0;
        goto exit;
    }
    

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto exit;
    }

    if (diagfac->count > 0) { diag_factor = diagfac->dval[0]; }
    if (gapfac->count > 0) { gap_factor = gapfac->dval[0]; }

    /* normal case: take the command line options at face value */
    exitcode = mymain(
            axes->count, binary->count,
            logbase->sval, logbase->count,
            images->count,
            local->count, outfiles->filename, outfiles->count,
            score->sval, score->count,
            diag_factor, gap_factor,
            testfile->filename, testfile->count, 
            infiles->filename, infiles->count,
            smat_out->filename, smat_out->count,
            smat_in->filename, smat_in->count
            );

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int mymain(int AXES, int BINARY, const char **logs, int nlogs, int IMAGES,
        int LOCAL, const char **outfiles, int outfile_cnt, 
        const char **SCORE_ARR, int score_cnt, 
        double diag_factor, double gap_factor,
        const char **testfile, int testfile_cnt, 
        const char **infiles, int ninfiles,
        const char **smat_out_files, int nsmat_out_files,
        const char **smat_in_files, int nsmat_in_files) {
  
    int NUM_INTERNAL_ANCHORS = 10000;   // Max this out!
    double _logbase;
    if (nlogs > 0) {
        _logbase = atof(logs[0]);
    }
    char file1[1024];
    char file2[1024];
    char outfile[1024];
    char SCORE[1024];
    char TESTFILE[1024];
    strcpy(file1, infiles[0]);
    strcpy(file2, infiles[1]);
    if (outfile_cnt) {
        strcpy(outfile, outfiles[0]);
    }
    if (score_cnt) {
        strcpy(SCORE, SCORE_ARR[0]);
    }
    else {
        strcpy(SCORE, "covariance");
    }
    if (testfile_cnt) {
        strcpy(TESTFILE, testfile[0]);
    }

    int mi_bins = 5;
    float init_penalty = 0.f;
    //char toError[300];

    if (strlen(SCORE) <= 1) {
        strcpy(SCORE, "covariance");
    }

    if (DEBUG) {
        std::cerr << "**********************************************\n";
        std::cerr << "LOCAL: " << LOCAL << "\n";
        std::cerr << "BINARY: " << BINARY << "\n";
        std::cerr << "IMAGES: " << IMAGES << "\n";
        std::cerr << "AXES: " << AXES << "\n";
        std::cerr << "SCORE: " << SCORE << "\n";
        std::cerr << "outfile: "<< outfile << "\n";
        std::cerr << "TESTFILE: " << TESTFILE << "\n";
        std::cerr << "file1: " << file1 << "\n";
        std::cerr << "file2: " << file2 << "\n";
        std::cerr << "**********************************************\n";
    }

    // ************************************************************
    // * READ IN FILES TO GET MAT 
    // ************************************************************
    LMat lmat1;
    LMat lmat2;
    MatF mat1;
    MatF mat2;
    MatF smat;
    DynProg dyn;

    if (AXES) {
        if (BINARY) {
            lmat1.set_from_binary(file1);
            lmat2.set_from_binary(file2);
            lmat1._mat->copy(mat1, 1);
            lmat2._mat->copy(mat2, 1);
        }
        else {
            lmat1.set_from_ascii(file1);
            lmat2.set_from_ascii(file2);
            lmat1._mat->copy(mat1, 1);
            lmat2._mat->copy(mat2, 1);
        }
    }
    else {
        //mat1.set_from_ascii(file1);  @TODO: write this guy
        //mat2.set_from_ascii(file2); 
    }

    // ************************************************************
    // * SCORE THE MATRICES
    // ************************************************************
    if (DEBUG) {
        std::cerr << "Scoring the mats!\n";
    }
    if (nsmat_in_files) {
        smat.set_from_binary(smat_in_files[0]);
        dyn._smat = &smat;
    }
    else {
        dyn.score(mat1, mat2, smat, SCORE, mi_bins);
        if (nlogs > 0) {
            smat.logarithm(_logbase);  
        }
        // SETTING THE SMAT TO BE IN TERMS OF PERCENT of AVERAGE!!!!!!
        smat *= 100/smat.avg();
    }
    if (nsmat_out_files) {
        printf("Writing binary smat to '%s'\n", smat_out_files[0]);
        smat.write(smat_out_files[0]);
        //smat.print(smat_out_files[0]);
        exit(0);
    }


    // ************************************************************
    // * PREPARE GAP PENALTY ARRAY
    // ************************************************************
   
    double average = smat.avg();
    int gp_length = smat.rows() + smat.cols();
    float max_gap_slope = 1000;
    //float max_gap_slope = 2;
    int num_b_steps = 1;
    float mstart = 0;

    int length_vecs = 0;
    float diff_ = max_gap_slope - mstart;
    float steps_ = diff_/15;
    for (float m = mstart; m <= max_gap_slope; m += steps_) {
       length_vecs++; 
    }

    // Write the file for plotting each guy:
    std::ofstream fh("linearopt");
    puts("WRITING TO LINEAROPT");
    fh << "gap_penalty_optimization" << "\n";
    fh << "linear gap penalty optimization (smat avg = " << smat.avg() << ")" << "\n";
    fh << "slope of linear gap penalty" << "\n";
    fh << "sum of the square residuals" << "\n";

    MatF tester;
    MatF tester_trans;
    VecF mpt;
    VecF npt;
    VecF mOut_tm;
    VecF nOut_tm;

    if (testfile_cnt) {
        tester.set_from_ascii(TESTFILE, 1);  // no headers on the files
        tester.transpose(tester_trans);
        mpt.set(tester_trans.cols(), tester_trans.pointer(0));
        npt.set(tester_trans.cols(), tester_trans.pointer(1));
    }

    for (float b = 0.f; b <= 1000.f; b=(1+b)*2) {
        // write the datalabel
        fh << "y intercept = " << b << "\n";
        VecD m_vals(length_vecs);
        VecD residual_vals(length_vecs);
        int vec_cnt = 0;
        for (float m = mstart; m <= max_gap_slope; m += steps_) {
            float slope = m;
            VecF gp_array;
            dyn.linear_less_before(m,b,gp_length,gp_array);

            printf("M: %f B: %f gplength %d\n", m,b,gp_length); //puts("GPARRAY: "); gp_array.print();

            // ************************************************************
            // * DYNAMIC PROGRAM
            // ************************************************************ 
            int minimize = 0;
            if (DEBUG) {
                std::cerr << "Dynamic Time Warping Score Matrix!\n";
            }
            if (!(strcmp(SCORE, "euclidean"))) {
                minimize = 1;
            }
            dyn.find_path(smat, gp_array, minimize, diag_factor, gap_factor, LOCAL, init_penalty);
            printf("DYNPROG SCORE: %f\n", dyn._bestScore);

            VecI mOut;
            VecI nOut;
            dyn.warp_map(mOut, nOut, minimize, NUM_INTERNAL_ANCHORS);

            //puts("MOUT:"); mOut.print(); puts("NOUT:"); nOut.print();

            if (testfile_cnt) {
                //puts("MOUT"); mOut.print();
                lmat1.tm_axis_vals(mOut, mOut_tm);
                lmat2.tm_axis_vals(nOut, nOut_tm);
                double sum_res = (double) (dyn.sum_sq_res_yeqx(mOut_tm, nOut_tm, mOut, nOut, mpt, npt));
                printf("sum res: %f\n", sum_res);
                //double old_abs_average = VecF::avg_abs_diff(npt, mpt);
                //double new_abs_average = VecF::avg_abs_diff(newnpt, mpt);
                //printf("ORIG AVG abs diff: %f\n", old_abs_average);
                //printf("NEW  AVG abs diff: %f\n", new_abs_average);

                // Set values for plotting:
                m_vals[vec_cnt] = slope;
                residual_vals[vec_cnt] = sum_res;
            }

            if (IMAGES) {
                PngIO wrt(1);
                char base_fn[1024];
                sprintf(base_fn, "i_%d", vec_cnt);
                char tb_fn[1024];
                strcpy(tb_fn, base_fn);
                strcat(tb_fn, "tb.png");
                //char *tb_fn = "tb.png";
                wrt.write(tb_fn, dyn._tb);
                char tbpath_fn[1024];
                strcpy(tbpath_fn, base_fn);
                strcat(tbpath_fn, "tbpath.png");
                wrt.write(tbpath_fn, dyn._tbpath);

                char asmat_fn[1024];
                strcpy(asmat_fn, base_fn);
                strcat(asmat_fn, "asmat.png");
                //wrt.write(asmat_fn, dyn._asmat);


                //strcpy(base_fn, "tb.png");
                //char *tbpath_fn = "tbpath.png";
                //char *tbscores_fn = "tbscores.png";
                //wrt.write(tbscores_fn, dyn._tbscores);
                //char *asmat_fn = "asmat.png";
                //wrt.write(asmat_fn, dyn._asmat);
                char *smat_fn = "smat.png";
                //wrt.write(smat_fn, *dyn._smat);
            }




            vec_cnt++;
        }
        m_vals.print(fh, 1); 
        residual_vals.print(fh, 1); 
    }
    fh.close();
    /*
    char silly[100];
    strcpy(silly, "png_");
    char tmpp[5];
    sprintf(tmpp, "%d", i);
    strcat(silly, tmpp); 
    strcat(silly, ".png");

    PngIO wrt(0);
    //wrt.write(silly, dyn._tbpath);
    wrt.write(silly, _scorepath);
    */

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

