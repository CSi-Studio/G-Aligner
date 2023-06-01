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

int mymain(int AXES, int BINARY, int IMAGES, int LOCAL, 
        double gap_init, double gap_extend,
        double factor_diag, double factor_gap,
        const char *smat_out_file, const char *smat_in_file,
        const char *timefile,
        const char *score, 
        const char **infiles, int infiles_cnt,
        const char **outfiles, int outfile_cnt 
        );


int main (int argc, char **argv) {
    // Create the argument structures:
    struct arg_lit  *axes_str     = arg_lit0("a", "axes", "file contains x and y axis labels");
    struct arg_lit  *binary_str   = arg_lit0("b", "binary", "binary data file");
    struct arg_lit  *help_str     = arg_lit0("h", "help", "prints this help and exits");
    struct arg_lit  *images_str   = arg_lit0("i", "images", "png images created");
    struct arg_lit  *local_str    = arg_lit0("l", "local", "local alignment");
    struct arg_file *outfiles_str = arg_file0("o", "outfile", "<outfile>", "output to file");
    struct arg_str  *score_str    = arg_str0("s", "score", "<scoretype>", "similarity score to compare vectors");
    struct arg_dbl  *gap_init_str = arg_dbl0(NULL, "gap_init", "<dbl>", "penalty for initiating a gap");
    struct arg_dbl  *gap_extend_str = arg_dbl0(NULL, "gap_extend", "<dbl>", "penalty for extending a gap");
    struct arg_dbl  *factor_diag_str  = arg_dbl0(NULL, "factor_diag", "2.0", "factor applied to diagonal steps");
    struct arg_dbl  *factor_gap_str = arg_dbl0(NULL, "factor_gap", "2.0", "factor applied to diagonal steps");
    struct arg_file *timefile_str = arg_file0("t", "timefile", "<times.txt>", "file with times to verify alignment");
    struct arg_file *infiles_str  = arg_filen(NULL, NULL, NULL, 2, 2, "files to align (first is template)");
    struct arg_file *smat_out_str = arg_file0(NULL, "smat_out", "<file>", "writes binary smat to file and exits");
    struct arg_file *smat_in_str  = arg_file0(NULL, "smat_in", "<file>", "uses binary smat file instead of calculating");
    struct arg_end  *end_str      = arg_end(20);

    void* argtable[] = {axes_str,binary_str,help_str,images_str,local_str,outfiles_str,score_str,gap_init_str,gap_extend_str,factor_diag_str,factor_gap_str,timefile_str,infiles_str,smat_out_str,smat_in_str,end_str};
    const char* progname = "obi-warp";
    int nerrors;
    int exitcode=0;
    char smat_in[1024];
    char smat_out[1024];
    char timefile[1024];



    // DECLARE default vars here:
    char score[20];
    double factor_diag = 2.f;
    double factor_gap = 1.f;
    double gap_init = 0.0;
    double gap_extend = 0.0;

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

    if (help_str->count > 0 || argc < 3)
    {
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
        printf("-s, --score = score function: covariance, product (dot product)\n"); 
        printf("     *pearsons_r, pearsons_r2, euclidean, mutual_info\n");
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
        arg_print_errors(stdout,end_str,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto exit;
    }

    // ******************************************************************
    // SET THE DEFAULTS:
    // ******************************************************************

    // score function:
    strcpy(score, "pearsons_r");  // DEFAULT SCORE
    if (score_str->count > 0) { strcpy(score, score_str->sval[0]); }

    // diag and gap factors:
    if (factor_diag_str->count > 0) { factor_diag = factor_diag_str->dval[0]; }
    if (factor_gap_str->count > 0) { factor_gap = factor_gap_str->dval[0]; }

    // gap penalty defaults (taken from nosize, residuals^2 summed results)
    if (!strcmp(score,"pearsons_r")) {
        gap_extend = 2.4;
        gap_init = 0.3;
    }
    else if (!strcmp(score,"covariance")) {
        gap_extend = 11.7;
        gap_init = 0.0;
    }
    else if (!strcmp(score,"euclidean")) {
        gap_extend = 1.8;
        gap_init = 0.9;
    }
    else if (!strcmp(score,"product")) {
        gap_extend = 7.8;
        gap_init = 0.0;
    }

    // gap values
    if (gap_init_str->count > 0) {
        gap_init = gap_init_str->dval[0];
    }
    if (gap_extend_str->count > 0) {
        gap_extend = gap_extend_str->dval[0];
    }

    strcpy(smat_in,"");
    if (smat_in_str->count > 0) {
        strcpy(smat_in, smat_in_str->filename[0]);
    }
    strcpy(smat_out,"");
    if (smat_out_str->count > 0) {
        strcpy(smat_out, smat_out_str->filename[0]);
    }

    strcpy(timefile, "");
    if (timefile_str->count > 0) {
        strcpy(timefile, timefile_str->filename[0]);
    }

    exitcode = mymain(
            axes_str->count, binary_str->count, 
            images_str->count, local_str->count, 
            gap_init, gap_extend,
            factor_diag, factor_gap,
            smat_in, smat_out, 
            timefile,
            score,
            infiles_str->filename, infiles_str->count,
            outfiles_str->filename, outfiles_str->count
            );

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}


int mymain(int AXES, int BINARY, int IMAGES, int LOCAL, 
        double gap_init, double gap_extend,
        double factor_diag, double factor_gap,
        const char *smat_out_file, const char *smat_in_file,
        const char *timefile,
        const char *score, 
        const char **infiles, int infiles_cnt,
        const char **outfiles, int outfile_cnt
        ) {
  
    bool have_smat_in_file = false;
    bool have_smat_out_file = false;
    bool have_timefile = false;
    if (strlen(smat_out_file) > 0) { have_smat_out_file = true; }
    if (strlen(smat_in_file) > 0) { have_smat_in_file = true; }
    if (strlen(timefile) > 0) { have_timefile = true; }

    int NUM_INTERNAL_ANCHORS = 10000;   // Max this out!
    char file1[1024];
    char file2[1024];
    char outfile[1024];
    strcpy(file1, infiles[0]);
    strcpy(file2, infiles[1]); 
    if (outfile_cnt) {
        strcpy(outfile, outfiles[0]);
    }


    int mi_bins = 5;
    float init_penalty = 0.f;
    //char toError[300];


    if (DEBUG) {
        std::cerr << "**********************************************\n";
        std::cerr << "LOCAL: " << LOCAL << "\n";
        std::cerr << "BINARY: " << BINARY << "\n";
        std::cerr << "IMAGES: " << IMAGES << "\n";
        std::cerr << "AXES: " << AXES << "\n";
        std::cerr << "SCORE: " << score << "\n";
        std::cerr << "outfile: "<< outfile << "\n";
        std::cerr << "TIMEFILE: " << timefile << "\n";
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
    if (have_smat_in_file) {
        smat.set_from_binary(smat_in_file);
        dyn._smat = &smat;
    }
    else {
        dyn.score(mat1, mat2, smat, score, mi_bins);
        // SETTING THE SMAT TO BE std normal
        smat -= smat.avg();
        double mean, stdev;
        smat._dat.sample_stats(mean, stdev);
        smat /= stdev;
        if (!strcmp(score,"euclidean")) {
            smat *= -1; // inverting euclidean
        }
    }
    if (have_smat_out_file) {
        printf("Writing binary smat to '%s'\n", smat_out_file);
        smat.write(smat_out_file);
        //smat.print(smat_out_files[0]);
        exit(0);
    }

    // ************************************************************
    // * PREPARE GAP PENALTY ARRAY
    // ************************************************************
   
    MatF tester;
    MatF tester_trans;
    VecF mpt;
    VecF npt;
    VecF mOut_tm;
    VecF nOut_tm;

    int gp_length = smat.rows() + smat.cols();

    VecF gp_array;
    dyn.linear_less_before(gap_extend,gap_init,gp_length,gp_array);

    // ************************************************************
    // * DYNAMIC PROGRAM
    // ************************************************************ 
    int minimize = 0;
    if (DEBUG) {
        std::cerr << "Dynamic Time Warping Score Matrix!\n";
    }
    dyn.find_path(smat, gp_array, minimize, factor_diag, factor_gap, LOCAL, init_penalty);


    char newscorename[100];
    strcpy(newscorename,score);
    if (!strcmp(newscorename,"pearsons_r")) {
        strcpy(newscorename, "pearsonsr");
    }

    // SETUP THE FILES TO PRINT:
    char resplotfilename[1024];
    strcpy(resplotfilename, timefile);
    strcat(resplotfilename, ".res.anchors_");
    strcat(resplotfilename, newscorename);
    char toplotresname[1024];
    strcpy(toplotresname, resplotfilename);
    strcat(toplotresname, ".toplot");

    std::ofstream fhres(toplotresname);
    printf("WRITING TO: %s\n", toplotresname);

    fhres << "XYData" << "\n"; // type of data
    fhres << resplotfilename << "\n";
    fhres << "number internal anchors optimization (normalized smat)\n";
    fhres << "num internal anchors\n";
    fhres << "SUM(res^2)\n";

    char avgresplotfilename[1024];
    strcpy(avgresplotfilename, timefile);
    strcat(avgresplotfilename, ".avgres.anchors_");
    strcat(avgresplotfilename, newscorename);
    char toplotavgresname[1024];
    strcpy(toplotavgresname, avgresplotfilename);
    strcat(toplotavgresname, ".toplot");

    std::ofstream fhavgres(toplotavgresname);
    printf("WRITING TO: %s\n", toplotavgresname);

    fhavgres << "XYData" << "\n"; // type of data
    fhavgres << avgresplotfilename << "\n";
    fhavgres << "number internal anchors optimization (normalized smat)\n";
    fhavgres << "num internal anchors\n";
    fhavgres << "AVG(res^2)\n";

    char absplotfilename[1024];
    strcpy(absplotfilename, timefile);
    strcat(absplotfilename, ".abs.anchors_");
    strcat(absplotfilename, newscorename);
    char toplotabsname[1024];
    strcpy(toplotabsname, absplotfilename);
    strcat(toplotabsname, ".toplot");

    std::ofstream fhabs(toplotabsname);
    printf("WRITING TO: %s\n", toplotabsname);

    fhabs << "XYData" << "\n"; // type of data
    fhabs << absplotfilename << "\n";
    fhabs << "number internal anchors optimization (normalized smat)\n";
    fhabs << "num internal anchors\n";
    fhabs << "avg abs delta (sec)\n";

    if (have_timefile) {
        tester.set_from_ascii(timefile, 1);  // no headers on the files
        tester.transpose(tester_trans);
        mpt.set(tester_trans.cols(), tester_trans.pointer(0));
        npt.set(tester_trans.cols(), tester_trans.pointer(1));
    }

    int mcnt = 0;
    int ncnt = 0;

    int num_its = 10000;

    int vec_cnt = 0;
    double *residual_vals_arr = new double[num_its];
    double *avgresidual_vals_arr = new double[num_its];
    double *avg_abs_diff_vals_arr = new double[num_its];
    int *xaxis_arr = new int[num_its];
 
    int num_internal_anchors_used;
    int prev_num_int_anchors = -1;
    for (int x = 0; x < num_its; ++x) {
        VecI mOut;
        VecI nOut;
        int num_int_anchors = x;
        num_internal_anchors_used = dyn.warp_map(mOut, nOut, minimize, num_int_anchors);
        if (num_internal_anchors_used == prev_num_int_anchors) {
            break;
        }
        xaxis_arr[x] = num_internal_anchors_used;

        if (have_timefile) {
            float ssr, asr, sad, aad;
            dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt, npt, ssr, asr, sad, aad);
            // Set values for plotting:
            residual_vals_arr[vec_cnt] = (double)ssr;
            avgresidual_vals_arr[vec_cnt] = (double)asr;
            avg_abs_diff_vals_arr[vec_cnt] = (double)aad;

        }
        vec_cnt++;
    }
    VecI xaxis(num_internal_anchors_used+1, xaxis_arr);
    VecD residual_vals(num_internal_anchors_used+1, residual_vals_arr);
    VecD avgresidual_vals(num_internal_anchors_used+1, avgresidual_vals_arr);
    VecD avg_abs_diff_vals(num_internal_anchors_used+1, avg_abs_diff_vals_arr);

    fhres << newscorename << "\n";
    xaxis.print(fhres,1);
    fhavgres << newscorename << "\n";
    xaxis.print(fhavgres,1);
    fhabs << newscorename << "\n";
    xaxis.print(fhabs,1);


    residual_vals.print(fhres, 1); 
    fhres.close();
    avgresidual_vals.print(fhavgres, 1); 
    fhavgres.close();
    avg_abs_diff_vals.print(fhabs, 1);
    fhabs.close();

return 0;
}

