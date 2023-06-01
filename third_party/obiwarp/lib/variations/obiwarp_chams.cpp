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

#ifndef min
	#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#endif
#ifndef max
	#define max(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#endif

void diff(VecF &first, VecF &sec, VecF &outFirstLessSecond);

void diff(VecF &first, VecF &sec, VecF &outFirstLessSecond) {
    first.copy(outFirstLessSecond);
    outFirstLessSecond -= sec;
}

// makes a vector for interpolation from the min,max of outline
void setup_interp(VecF &outline, VecF &out);
void setup_interp(VecF &outline, VecF &out) {
    float outline_min, outline_max;
    outline.min_max(outline_min, outline_max);
    // make x vector of all the m points in between:
    int start_m = (int)outline_min;
    int end_m = 1 + (int)outline_max;
    int num_pts = end_m - start_m + 1;
    VecF out_tmp(num_pts);
    int _cnt = 0;
    for (int i = start_m; i <= end_m; ++i) {
        out_tmp[_cnt] = i; 
        _cnt++;
    }
    out.take(out_tmp);
}


int mymain(int AXES, int BINARY, int IMAGES, int LOCAL, 
        double gap_init, double gap_extend,
        double factor_diag, double factor_gap,
        const char *smat_out_file, const char *smat_in_file,
        const char *timefile,
        const char *chamsfile,
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
    struct arg_file *chamsfile_str = arg_file0("c", "chamsfile", "<file.mychams>", "my chams alignment file");
    struct arg_file *infiles_str  = arg_filen(NULL, NULL, NULL, 2, 2, "files to align (first is template)");
    struct arg_file *smat_out_str = arg_file0(NULL, "smat_out", "<file>", "writes binary smat to file and exits");
    struct arg_file *smat_in_str  = arg_file0(NULL, "smat_in", "<file>", "uses binary smat file instead of calculating");
    struct arg_end  *end_str      = arg_end(20);

    void* argtable[] = {axes_str,binary_str,help_str,images_str,local_str,outfiles_str,score_str,gap_init_str,gap_extend_str,factor_diag_str,factor_gap_str,timefile_str, chamsfile_str, infiles_str,smat_out_str,smat_in_str,end_str};
    const char* progname = "obi-warp";
    int nerrors;
    int exitcode=0;
    char smat_in[1024];
    char smat_out[1024];
    char timefile[1024];
    char chamsfile[1024];



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

    strcpy(chamsfile, "");
    if (chamsfile_str->count > 0) {
        strcpy(chamsfile, chamsfile_str->filename[0]);
    }

    exitcode = mymain(
            axes_str->count, binary_str->count, 
            images_str->count, local_str->count, 
            gap_init, gap_extend,
            factor_diag, factor_gap,
            smat_in, smat_out, 
            timefile,
            chamsfile,
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
        const char *chamsfile,
        const char *score, 
        const char **infiles, int infiles_cnt,
        const char **outfiles, int outfile_cnt
        ) {
  
    bool have_smat_in_file = false;
    bool have_smat_out_file = false;
    bool have_timefile = false;
    bool have_chamsfile = false;
    if (strlen(smat_out_file) > 0) { have_smat_out_file = true; }
    if (strlen(smat_in_file) > 0) { have_smat_in_file = true; }
    if (strlen(timefile) > 0) { have_timefile = true; }
    if (strlen(chamsfile) > 0) { have_chamsfile = true; }

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
        std::cerr << "TIMEFILE: " << chamsfile << "\n";
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

    double average = smat.avg();
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

    VecI mOut;
    VecI nOut;
    dyn.warp_map(mOut, nOut, minimize, NUM_INTERNAL_ANCHORS);

    VecF mOutFt;
    VecF nOutFt;

    lmat1.tm_axis_vals(mOut, mOutFt);
    lmat2.tm_axis_vals(nOut, nOutFt);

    if (have_chamsfile && have_timefile) {
        // LOAD TIMEFILE STUFF:
        tester.set_from_ascii(timefile, 1);  // no headers on the files
        tester.transpose(tester_trans);
        mpt.set(tester_trans.cols(), tester_trans.pointer(0));
        npt.set(tester_trans.cols(), tester_trans.pointer(1));

        // CHAMSFILE STUFF:
        MatF chams;
        chams.set_from_ascii(chamsfile, 1);
        int num_rows;
        VecF *vecs = new VecF[chams.rows()];
        chams.row_vecs(num_rows, vecs);

        // LOAD the chams file:
        VecI mCham, nCham;
        VecF mChamFt, nChamFt;
        mChamFt.set(vecs[0]);
        nChamFt.set(vecs[1]);
        vecs[2].to_i(mCham);
        vecs[3].to_i(nCham);
        float minMChamFt, maxMChamFt;
        mChamFt.min_max(minMChamFt, maxMChamFt);

        // Remove all mpt vals that are < or > range of chams
        VecF mpt_short;
        VecF npt_short;
        mpt.copy(mpt_short);
        npt.copy(npt_short);
        int *sh_indices_arr = new int[mpt_short.size()];
        int num_too_short = 0;
        for (int i = 0; i < mpt_short.size(); ++i) {
            if ( (mpt_short[i] < minMChamFt) || (mpt_short[i] > maxMChamFt)) {
                sh_indices_arr[num_too_short] = i;
                num_too_short++;
            }
        }
        VecI sh_indices(num_too_short, sh_indices_arr);
        sh_indices.sort();
        // Need to remove indices in reverse order so they have efficacy:
        for (int i = sh_indices.size() - 1; i >=0; --i) {
            mpt_short.remove(sh_indices[i]);
            npt_short.remove(sh_indices[i]);
        }
        //printf("BEFORE: %d %d AFTER: %d %d\n", mpt.size(), npt.size(), mpt_short.size(), npt_short.size()); mpt_short.print(); npt_short.print(); exit(1);


        // *************************************************
        // the m plot will go on the X axis
        // n on the Y axis
        // *************************************************
        // SETUP PLOT FILE:
        char fnplot[1024];
        char fnplot_toplot[1024];
        strcpy(fnplot, chamsfile);
        char *ptr = strrchr(fnplot, '/');
        strcpy(fnplot, ++ptr);
        // remove the mychams at the end
        ptr = strstr(fnplot, ".mychams");
        *ptr = '\0';
        strcpy(fnplot_toplot, fnplot);
        strcat(fnplot_toplot, ".toplot");
        std::ofstream fhplot(fnplot_toplot);
        fhplot << "XYData" << "\n";
        fhplot << fnplot << "\n";
        fhplot << "chams comparison" << "\n";
        fhplot << "run1 (sec)" << "\n";
        fhplot << "run2 (sec)" << "\n";

        // Sqaure residuals
        char fn_sr_indplot[1024];
        char fn_sr_indplot_toplot[1024];
        strcpy(fn_sr_indplot, fnplot);
        strcat(fn_sr_indplot, ".sr");
        strcpy(fn_sr_indplot_toplot, fn_sr_indplot);
        strcat(fn_sr_indplot_toplot, ".toplot");
        std::ofstream fh_sr_plot(fn_sr_indplot_toplot);
        fh_sr_plot << "YData" << "\n";
        fh_sr_plot << fn_sr_indplot << "\n";
        fh_sr_plot << "chams comparison" << "\n";
        fh_sr_plot << "MS runs" << "\n";
        fh_sr_plot << "res^2 (sec^2)" << "\n";

        // Absolute differences
        char fn_ad_indplot[1024];
        char fn_ad_indplot_toplot[1024];
        strcpy(fn_ad_indplot, fnplot);
        strcat(fn_ad_indplot, ".ad");
        strcpy(fn_ad_indplot_toplot, fn_ad_indplot);
        strcat(fn_ad_indplot_toplot, ".toplot");
        std::ofstream fh_ad_plot(fn_ad_indplot_toplot);
        fh_ad_plot << "YData" << "\n";
        fh_ad_plot << fn_ad_indplot << "\n";
        fh_ad_plot << "chams comparison" << "\n";
        fh_ad_plot << "MS runs" << "\n";
        fh_ad_plot << "avg abs diff (sec)" << "\n";

        // PRINT THE HEADER FOR THE OUTPUT
        // ORDER: chfemax, chfemin, linmax, linmin
        printf("\tNO_Align\t\t\t\tNO_Align_short\t\t\t\t%s_M%.1fB%.1f\t\t\t\t%s_M%.1fB%.1f_short\t\t\t\tchamsMaxChfe\t\t\t\tchamsMinChfe\t\t\t\tchamsMaxLinear\t\t\t\tchamsMinLinear\n", score, gap_extend, gap_init, score, gap_extend, gap_init);
        int num_datasets = 8;
        for (int i = 0; i < num_datasets; ++i) {
            printf("\tssr\tasr\tsad\taad");
        }
        printf("\tNUM_REMOVED\n");
        printf("%s\t", fnplot);
        // Path accuracy
        float ssr, asr, sad, aad;

        // Path accuracy with NO alignment
        // Create a linear warp path to test
        float min_ft = min(mOut[0], nOut[0]);
        float max_ft = max(mOut.last(), nOut.last());
        VecF lin_ft(2);
        lin_ft[0] = min_ft; lin_ft[1] = max_ft;
        dyn.path_accuracy(lin_ft, lin_ft, mpt, npt, ssr, asr, sad, aad, 1);
        printf("%f\t%f\t%f\t%f\t", ssr, asr, sad, aad);
        
        // Path accuracy with NO alignment (short)
        dyn.path_accuracy(lin_ft, lin_ft, mpt_short, npt_short, ssr, asr, sad, aad, 1);
        printf("%f\t%f\t%f\t%f\t", ssr, asr, sad, aad);

        // NO Align: DETAILS
        VecF no_align_sr;
        VecF no_align_ad;
        dyn.path_accuracy_details(lin_ft, lin_ft, mpt, npt, no_align_sr, no_align_ad, 1);
        fh_ad_plot << fn_ad_indplot << "_no_align\n";
        fh_sr_plot << fn_sr_indplot << "_no_align\n";
        no_align_sr.print(fh_sr_plot, 1);
        no_align_ad.print(fh_ad_plot, 1);

        // NO Align: DETAILS (short)
        VecF no_align_sr_short;
        VecF no_align_ad_short;
        dyn.path_accuracy_details(lin_ft, lin_ft, mpt_short, npt_short, no_align_sr_short, no_align_ad_short, 1);
        fh_ad_plot << fn_ad_indplot << "_no_align_short\n";
        fh_sr_plot << fn_sr_indplot << "_no_align_short\n";
        no_align_sr_short.print(fh_sr_plot, 1);
        no_align_ad_short.print(fh_ad_plot, 1);


        // Complete path accuracy:
        dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt, npt, ssr, asr, sad, aad);
        printf("%f\t%f\t%f\t%f\t", ssr, asr, sad, aad);

        // Short path accuracy
        dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt_short, npt_short, ssr, asr, sad, aad);
        printf("%f\t%f\t%f\t%f", ssr, asr, sad, aad);

        // MINE: DETAILS
        VecF mine_sr_short;
        VecF mine_ad_short;
        VecF mine_sr;
        VecF mine_ad;
        
        // NORMAL
        dyn.path_accuracy_details(mOutFt, nOutFt, mpt, npt, mine_sr, mine_ad, 0);
        fh_ad_plot << fn_ad_indplot << "_mine\n";
        fh_sr_plot << fn_sr_indplot << "_mine\n";
        mine_sr.print(fh_sr_plot, 1);
        mine_ad.print(fh_ad_plot, 1);

        // Short
        dyn.path_accuracy_details(mOutFt, nOutFt, mpt_short, npt_short, mine_sr_short, mine_ad_short, 0);
        fh_ad_plot << fn_ad_indplot << "_mine_short\n";
        fh_sr_plot << fn_sr_indplot << "_mine_short\n";
        mine_sr_short.print(fh_sr_plot, 1);
        mine_ad_short.print(fh_ad_plot, 1);


        // PLOT THE CHAMS PATH:
        fhplot << "chams" << "\n";
        mChamFt.print(fhplot, 1);
        //VecF nChamFt_minus_mChamFt;
        //diff(nChamFt, mChamFt, nChamFt_minus_mChamFt);
        nChamFt.print(fhplot, 1);

        // PLOT THE INTERPOLATED CHAMS PATHS:
        VecI mOutChams;
        VecI nOutChams;
        for (int linear = 0; linear < 2; ++linear) {
            for (minimize = 0; minimize < 2; ++minimize) {
                dyn.warp_map(mCham, nCham, vecs[4], mOutChams, nOutChams, minimize, NUM_INTERNAL_ANCHORS);

                VecF m_warp_times(mOutChams.length());
                VecF n_warp_times(nOutChams.length());
                for (int i = 0; i < mOutChams.length(); ++i) {
                    m_warp_times[i] = mChamFt[mCham.index(mOutChams[i])];
                    n_warp_times[i] = nChamFt[nCham.index(nOutChams[i])];
                }
                //puts("mwarp"); m_warp_times.print(); puts("nwarp"); n_warp_times.print(); exit(1);

                float ssr, asr, sad, aad;
                VecF chams_sr_plot;
                VecF chams_ad_plot;
                dyn.path_accuracy_details(m_warp_times, n_warp_times, mpt_short, npt_short, chams_sr_plot, chams_ad_plot, linear);
                dyn.path_accuracy(m_warp_times, n_warp_times, mpt_short, npt_short, ssr, asr, sad, aad, linear);
                printf("\t%f\t%f\t%f\t%f", ssr, asr, sad, aad);

                if (minimize) {  // minimize 
                    if (linear) {
                        fhplot << "chams min linear" << "\n";
                        fh_ad_plot << fn_ad_indplot << "_min_lin\n";
                        fh_sr_plot << fn_sr_indplot << "_min_lin\n";
                    }
                    else {
                        fhplot << "chams min chfe" << "\n";
                        fh_ad_plot << fn_ad_indplot << "_min_interp\n";
                        fh_sr_plot << fn_sr_indplot << "_min_interp\n";
                    }
                }
                else {  // maximize
                    if (linear) {
                        fhplot << "chams max linear" << "\n";
                        fh_ad_plot << fn_ad_indplot << "_max_lin\n";
                        fh_sr_plot << fn_sr_indplot << "_max_lin\n";
                    }
                    else {
                        fhplot << "chams max chfe" << "\n";
                        fh_ad_plot << fn_ad_indplot << "_max_interp\n";
                        fh_sr_plot << fn_sr_indplot << "_max_interp\n";
                    }
                }
                // see interpolation
                VecF chams_m_interp;
                setup_interp(m_warp_times, chams_m_interp);
                VecF n_out_interp;
                if (linear) {
                    VecF::linear_interp(m_warp_times, n_warp_times, chams_m_interp, n_out_interp);
                }
                else {
                    VecF::chfe(m_warp_times, n_warp_times, chams_m_interp, n_out_interp);
                }
                //VecF n_out_interp_minus_m;
                //diff(n_out_interp, m_interp, n_out_interp_minus_m); 
                chams_m_interp.print(fhplot, 1);
                n_out_interp.print(fhplot, 1);
                chams_ad_plot.print(fh_ad_plot, 1);
                chams_sr_plot.print(fh_sr_plot, 1);
            }
        }
        // FINISH THE OUTPUT
        printf("\t%d", num_too_short);
        puts("");

        // PLOT MY PATH:
        fhplot << "path (default params)" << "\n";
        mOutFt.print(fhplot, 1);
        //VecF nOutFt_minus_mOutFt;
        //diff(nOutFt, mOutFt, nOutFt_minus_mOutFt);
        nOutFt.print(fhplot, 1);

        // SETUP THE INTERPOLANT POINTS
        VecF my_m_interp;
        setup_interp(mOutFt, my_m_interp);
        
        // PLOT MY INTERPOLANT:
        fhplot << "my interp" << "\n";
        VecF my_n_out_interp;
        VecF::chfe(mOutFt, nOutFt, my_m_interp, my_n_out_interp);
        my_m_interp.print(fhplot, 1);
        my_n_out_interp.print(fhplot, 1);

        // PLOT THE TIME PTS:
        fhplot << "time standards" << "\n"; 
        mpt.print(fhplot, 1);
        //VecF npt_minus_mpt;
        //diff(npt, mpt, npt_minus_mpt);
        npt.print(fhplot, 1);
        
        fhplot.close();
        fh_sr_plot.close();
        fh_ad_plot.close();
    }

    // Warp the second lmat run!
    if (AXES) {
        VecF nOutF;
        VecF mOutF;
        lmat1.tm_axis_vals(mOut, mOutF);
        lmat2.tm_axis_vals(nOut, nOutF); //
        lmat2.warp_tm(nOutF, mOutF); 
    }
    else {
        // or warp the mat itself!
        // @TODO: write the warping of mat itself!
    }

    char outfilename[1024 + 7];
    strcpy(outfilename, file2);
    strcat(outfilename, ".warped");
    if (BINARY) {
        //lmat2.print(outfilename);
        lmat2.write(outfilename);
    }
    else {
        lmat2.print(outfilename);
    }
    


    if (IMAGES) {
        PngIO wrt(1);
        char base_fn[1024];
        strcpy(base_fn, "align_");
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

