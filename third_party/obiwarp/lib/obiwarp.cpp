// STDLIB:
#include <cstdio>
#include <iostream>
#include <fstream>
#include "string.h"
#include <cstdlib>

// MINE
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "dynprog.h"
#include "pngio.h"
#include "cmdparser.h"

/********************************************/
char * VERSION = (char *)"0.9.4";
/********************************************/

#define DEBUG (0)

bool format_is_labelless(const char *format);

int main (int argc, char **argv) {
    // NOTE: use outfile as indicator if option passed in as opts.outfile!
    // because we can set opts.outfile to NULL and other routines will
    // automatically write to stdout!
    bool outfile = 0;
    bool outfile_is_stdout = 0;

    CmdParser opts(argc, argv, VERSION);

    if (opts.outfile != NULL) {
        outfile = 1;
        if (!strcmp(opts.outfile, "STDOUT")) {
            outfile_is_stdout = 1;
            opts.outfile = NULL;
        }
    }

    char file1[1024];
    char file2[1024];
    strcpy(file1, opts.infiles[0]);
    strcpy(file2, opts.infiles[1]); 

    // ************************************************************
    // * READ IN FILES TO GET MAT 
    // ************************************************************
    LMat lmat1;
    LMat lmat2;
    MatF smat;
    DynProg dyn;

    if (!strcmp(opts.format, "mat")) {
        lmat1.set_from_binary_mat(file1);
        lmat2.set_from_binary_mat(file2);
    }
    else if (!strcmp(opts.format, "mata")) {
        lmat1.set_from_ascii_mat(file1);
        lmat2.set_from_ascii_mat(file2);
    }
    else if (!strcmp(opts.format, "lmat")) {
        lmat1.set_from_binary(file1);
        lmat2.set_from_binary(file2);
    }
    else if (!strcmp(opts.format, "lmata")) {
        lmat1.set_from_ascii(file1);
        lmat2.set_from_ascii(file2);
    }
    ////puts("LMAT1 AND LMAT2"); lmat1.print(); lmat2.print();

    // ************************************************************
    // * SCORE THE MATRICES
    // ************************************************************
    if (DEBUG) {
        std::cerr << "Scoring the mats!\n";
    }
    if (opts.smat_in != NULL) {
        smat.set_from_binary(opts.smat_in);
        dyn._smat = &smat;
    }
    else {
        dyn.score(*(lmat1.mat()), *(lmat2.mat()), smat, opts.score);
        // SETTING THE SMAT TO BE std normal
        if (!opts.nostdnrm) {
            if (!smat.all_equal()) { 
                smat.std_normal();
            }
        }
        if (!strcmp(opts.score,"euc")) {
            smat *= -1; // inverting euclidean
        }
    }
    if (opts.smat_out != NULL) {
        std::cerr << "Writing binary smat to '" << opts.smat_out << "'\n";
        smat.write(opts.smat_out);
        //smat.print(smat_out_files[0]);
        exit(0);
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
    dyn.linear_less_before(opts.gap_extend,opts.gap_init,gp_length,gp_array);

    // ************************************************************
    // * DYNAMIC PROGRAM
    // ************************************************************ 
    int minimize = 0;
    if (DEBUG) {
        std::cerr << "Dynamic Time Warping Score Matrix!\n";
    }
    dyn.find_path(smat, gp_array, minimize, opts.factor_diag, opts.factor_gap, opts.local, opts.init_penalty);

    VecI mOut;
    VecI nOut;
    dyn.warp_map(mOut, nOut, opts.response, minimize);
    //puts("mOUT"); mOut.print(); nOut.print();

    // Major output unless its the only case where we don't need warped time
    // values
    if (!(outfile_is_stdout && format_is_labelless(opts.format))) {
        // MAJOR OUTPUT:
        VecF nOutF;
        VecF mOutF;
        lmat1.tm_axis_vals(mOut, mOutF);
        lmat2.tm_axis_vals(nOut, nOutF); //
        lmat2.warp_tm(nOutF, mOutF); 
        lmat2.tm()->print(1);
    }

    // No labels on matrix and we have an outfile to produce
    // Needs to be after MAJOR OUTPUT since it warps the data!
    if (format_is_labelless(opts.format) && outfile) {
        // @TODO: implement data warping here
    }

    // All subroutines below should write to the specified file
    // if the file == NULL then they should write to stdout!
    // opts.outfile is set to NULL if "STDOUT" is specified!
    if (outfile) {
        if (!strcmp(opts.format, "mat")) {
            lmat2.mat()->write(opts.outfile);
        }
        else if (!strcmp(opts.format, "mata")) {
            lmat2.mat()->print(opts.outfile);
        }
        else if (!strcmp(opts.format, "lmat")) {
            lmat2.write(opts.outfile);
        }
        else if (!strcmp(opts.format, "lmata")) {
            lmat2.print(opts.outfile);
        }
        else {
            std::cerr << "Can't output to" << opts.format << "format (yet)\n";
            exit(0);
        }
    }

    // After all other output to stdout
    if (opts.timefile != NULL) {
        time_tester.set_from_ascii(opts.timefile, 1);  // no headers on the files
        time_tester.transpose(time_tester_trans);
        mpt.set(time_tester_trans.cols(), time_tester_trans.pointer(0));
        npt.set(time_tester_trans.cols(), time_tester_trans.pointer(1));
        float ssr, asr, sad, aad;
        dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt, npt, ssr, asr, sad, aad);
        printf("%f %f %f %f\n", ssr, asr, sad, aad);
    }


    if (opts.images) {
        PngIO wrt(1);
        char base_fn[1024];
        strcpy(base_fn, "obi-warp_");
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
        char *smat_fn = (char *)"smat.png";
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

bool format_is_labelless(const char *format) {
    if (!strcmp(format,"mat") || !strcmp(format,"mata")) {
        return 1;
    }
    else {
        return 0;
    }
}
